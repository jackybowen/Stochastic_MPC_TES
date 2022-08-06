module Controller
include("opt_engine.jl")
using .Opt_engine

using CSV
using DataFrames
using Statistics
using JSON
# using Dates
# using MAT

# clearconsole()
function run(GUI_input::Dict)

    ###### read inputs ######

    for (key, value) in GUI_input
        @eval $(Symbol(key)) = $value
    end
	###### Inputs ######
	Delta_T = 60;	 	# prediction resolution in min,
	dT = Delta_T/60;    # Delta_T in hour
	dT_m = 5;           # data resolutin in min,
	H = Int(24/dT);
	N = 14;             # Max length of simulation in days
	Ndata = N*24*60/dT_m# Max number of baseline data points
	# Tr = 24;     		# Tzon Setpoint in C
	solverflag = 1; 	# 1 for the 1st-order RC model
	uncertoat_flag = 0; # 0 for determinstic case of Toa, 1 for stochastic case of Toa
	uncertout_flag = 0; # 0 for determinstic case of Tout, 1 for stochastic case of Tout
	optcost_flag = 1;	# enable price profile in the cost function

	Tbd = 2; 			# bandwidth of comfort bounds around Tsetpoint

	numzones = 25		# Number of zones in building
	numahu = 4			# Number of AHUs in building
	# AHU serve zone index: AHU1: 1~17 AHU3:18~23 AHU2:24, AHU4:25
	ahuzonelist = [1:17,24,18:23,25]
	mzonelist = [1:23;]
	szonelist = [24,25]
	mahulist = [1,3]
	sahulist = [2,4]

	StartTime = 3300
	CurrentTime = Int64(Time[end])
	if CurrentTime <= StartTime
		t = 60;
	else
		t = Int(ceil.((CurrentTime - StartTime)/60))+60;
	end
	println("Current time is Hour ",Int(t/60))

##  Receive measurement at t
	zonelist  = ["VAV-102", "VAV-118", "VAV-119","VAV-120","VAV-123A","VAV-123B","VAV-127A","VAV-127B","VAV-129","VAV-131","VAV-133","VAV-136","VAV-142","VAV-143","VAV-150","VAV-CORRIDOR","VAV-RESTROOM","VAV-104","VAV-105","VAV-107","VAV-108","VAV-112","VAV-116","AHU-002","AHU-004"]
	tempsetpoints = ["ZONE-$z:Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)" for z in zonelist]
	tdiss = ["AHU-00$f SUPPLY EQUIPMENT OUTLET NODE:System Node Temperature [C](TimeStep)" for f in 1:numahu]
	T_oa = mean(eval(OutdoorAirTemperature));
	Tinit = zeros(numzones);
	mflow_init = zeros(numzones);
	for z = 1:numzones
		Tinit[z] = mean(eval(Symbol("Tzon$z")))
		mflow_init[z] = mean(eval(Symbol("Mflow$z")))
	end
##  Load data from baseline csv files
	t_i = Int(t/dT_m);
	Tr = zeros(numzones,H)
	T_oa_p = zeros(H)
	bsl_tempdischarge = zeros(numahu,H)
	for i=1:H
		thrz = t_i+Int((i-1)*60/dT_m)+1:t_i+Int(i*60/dT_m)
		if thrz[end] > Ndata
			thrz = Int.(mod1.(thrz, Ndata)) # Cycle indices when exceeding the max length of data
		end
		for z=1:numzones
			Tr[z,i] = mean(CSV.read(string("profile/baseline_setpoint.csv"),DataFrame)[!,Symbol(tempsetpoints[z])][thrz]);
		end
		for f = 1:numahu
			bsl_tempdischarge[f,i] = mean(CSV.read(string("profile/baseline_tempdischarge.csv"),DataFrame)[!,Symbol(tdiss[f])][thrz]);
		end
		T_oa_p[i] = mean(CSV.read(string("profile/baseline_oat.csv"),DataFrame)[!,Symbol("Environment:Site Outdoor Air Drybulb Temperature [C](TimeStep)")][thrz]);
	end
	Tzmin = Tr.+Tbd;
	Tzmax = Tr.-Tbd;
	println("Hourly OutdoorAirTemperatures for next 24 hours: ", T_oa_p)

	# T_oa_p = repeat(Toa_p,inner=(Int(60/Delta_T),1));
	price = CSV.read(string("profile/daily_prices.csv"), DataFrame)[!,Symbol("price")][1:end];# Minutely price for one day
	price[721:1080] = price[721:1080].*2; # Double the price during occupied period(12pm to 18pm)
	min_idx = repeat(1:Int(24*60), N+1);

##  Read  Model parameters from JSON files
	zone_coeffs = zeros(numzones,4)
	fan_coeffs = zeros(numahu,4)
	chiller_coeffs = zeros(3)
	# Pm = zeros(numzones)
	for z = 1:numzones
		coeff_dict = JSON.parsefile("profile/model_params/coeff_$(zonelist[z]).json")
		zone_coeffs[z,1] = coeff_dict["temp_pre"]
		zone_coeffs[z,2] = coeff_dict["tout"]
		zone_coeffs[z,3] = coeff_dict["tempset"]
		zone_coeffs[z,4] = coeff_dict["intercept"]
	end
	for f = 1:numahu
		coeff_dict = JSON.parsefile("profile/model_params/coeff_ahu_$f.json")
		fan_coeffs[f,1] = coeff_dict["flow"]
		fan_coeffs[f,2] = coeff_dict["flow2"]
		fan_coeffs[f,3] = coeff_dict["flow3"]
		fan_coeffs[f,4] = coeff_dict["intercept"]
	end
	coeff_dict = JSON.parsefile("profile/model_params/coeff_chiller.json")
	chiller_coeffs[1] = coeff_dict["load"]
	chiller_coeffs[2] = coeff_dict["load2"]
	chiller_coeffs[3] = coeff_dict["intercept"]
	# # Tinit = repeat(Tinit, Int(numzones/10));
	# global u_t = zeros(numzones)
	# global Opcost = zeros(1,Ntime)
	global Mf_p = zeros(numzones)
	global Tset_p = zeros(numzones)
	# global cputime = 0;

	global Optimization_input = Dict("numzones" => numzones,
									"numahu"    =>numahu,
									"ahuzonelist"=>ahuzonelist,
									"mzonelist" => mzonelist,
									"szonelist" => szonelist,
									"mahulist"  => mahulist,
									"sahulist"  => sahulist,
									"coeffs"    =>zone_coeffs,
									"fan_params"=>fan_coeffs,
									"chiller_params"=>chiller_coeffs,
									"price" 	=> price,
									"min_idx"	=>min_idx,
									# "hour_idx"=>hour_idx,
									"dT" 	=> dT,
									"t_0"	=> t,
									"H"		=> H,
									"T_init"=>Tinit,
									"m_init"=>mflow_init,
									"T_oa"  => T_oa,
									"T_oa_p"=> T_oa_p,
									"Tzmin" => Tzmin,
									"Tzmax" => Tzmax,
									# "ahudischargetemp_min"=>0,
									# "ahudischargetemp_max"=>30,
									"dischargetemp" =>bsl_tempdischarge,
									"damper_min" 	=> 0,
									"damper_max" 	=> 1,
									"zoneflow_min" 	=> 0,
									"zoneflow_max" 	=> 2,
									"solverflag"	=>solverflag,
									"uncertoat_flag"=>uncertoat_flag,
									"uncertout_flag"=>uncertout_flag,
									"optcost_flag"  =>optcost_flag);
	# "sto_model_flag"=>sto_model_flag,
	# "numrl" 		=>numrl);

	# ###### call optimization ######
	global sol = Opt_engine.run(Optimization_input)
	# println(cputime)
	for z = 1:numzones
		Mf_p[z] = sol[!,r"mflow"][1,z];
		Tset_p[z] = sol[!,r"tzon"][1,z];
		# Mf_p_store[z,Int(ceil(t/60/dT_m))] = Mf_p[z];
	end
	# end

	println("Next massflow target:  \n\t",Mf_p)
	println("Output Temperature setpoint: \n\t",Tset_p)
	# uncontrolled_power = CSV.read(string("rest_loadpower.csv"), DataFrame)[!,Symbol("Hourlypower")][Int(ceil(t/60))];
	# result = Dict("demand_target"=>sum(Mf_p[z] for z=1:numzones)+uncontrolled_power)

	CSV.write("saved/sol_T$t.csv", sol, append=false)
	result = Dict("Tset$z"=> Tset_p[z] for z=1:numzones)
	return result

end
end
