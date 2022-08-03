module Controller
include("opt_engine.jl")
using .Opt_engine

using CSV
using DataFrames
using Statistics
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
	N = 31;             # Max length of simulation in days
	# Tr = 24;     		# Tzon Setpoint in C
	solverflag = 1; 	# 1 for the 1st-order RC model
	uncertoat_flag = 0; # 0 for determinstic case of Toa, 1 for stochastic case of Toa
	uncertout_flag = 0; # 0 for determinstic case of Tout, 1 for stochastic case of Tout
	optcost_flag = 1;	# enable price profile in the cost function

	Tbd = 2; 			# bandwidth of T bounds around Tsetpoint

	numzones = 25		# Number of zones in building
	numahu = 4			# Number of AHUs in building
	numcoeffs = 4 		# Number of coefficients of zone-temp model
	# AHU serve zone index: AHu1: 1~17 AHU3:18~23 AHU2:24, AHU4:25
	ahuzonelist = [1:17,24,18:23,25]

	StartTime = 3300
	CurrentTime = Int64(Time[end])
	println("Current time is ",CurrentTime)
	if CurrentTime <= StartTime
		t = 60;
	else
		t = Int(ceil.((CurrentTime - StartTime)/60))+60;
	end

	zonelist  = ["VAV-102", "VAV-118", "VAV-119","VAV-120","VAV-123A","VAV-123B","VAV-127A","VAV-127B","VAV-129","VAV-131","VAV-133","VAV-136","VAV-142","VAV-143","VAV-150","VAV-CORRIDOR","VAV-RESTROOM","VAV-104","VAV-105","VAV-107","VAV-108","VAV-112","VAV-116","AHU-002","AHU-004"]
	tempsetpoints = ["ZONE-$z:Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)" for z in zonelist]
	T_oa = Float64(OutdoorAirTemperature[end]);
	Tinit = zeros(numzones);
	for z = 1:numzones
		Tinit[z] = mean(eval(Symbol("Tzon$z")))
	end

	t_i = Int(t/dT_m);
	Tr = zeros(numzones,H)
	T_oa_p = zeros(H)
	Tdischarge = zeros(numahu,H)
	for i=1:H
		for z=1:numzones
			Tr[z,i] = CSV.read(string("profile/baseline_setpoint.csv"),DataFrame)[!,Symbol(tempsetpoints[z])][t_i+Int((i-1)*60/dT_m)];
		end
		for f = 1:numahu
			Tdischarge[f,i] = CSV.read(string("profile/baseline_tempdischarge.csv"),DataFrame)[!,Symbol("AHU$f")][t_i+Int((i-1)*60/dT_m)];
		end
		T_oa_p[i] = CSV.read(string("profile/baseline_oat.csv"),DataFrame)[!,Symbol("Environment:Site Outdoor Air Drybulb Temperature [C](TimeStep)")][t_i+Int((i-1)*60/dT_m)];
	end
	Tzmin = Tr.+Tbd;
	Tzmax = Tr.-Tbd;

	# T_oa_p = repeat(Toa_p,inner=(Int(60/Delta_T),1));
	price = CSV.read(string("profile/daily_prices.csv"), DataFrame)[!,Symbol("price")][1:end];# Minutely price for one day
	price[721:1080] = price[721:1080].*2; # Double the price during occupied period(12pm to 18pm)
	min_idx = repeat(1:Int(24*60), N+1);

	Coeffs_file1 = string("profile/zonetemp_params.csv");
	Coeffs_file2 = string("profile/fan_sepa_params.csv");
	Coeffs_file3 = string("profile/chiller_sepa_params.csv");
	coeffs = zeros(numzones,numcoeffs)
	fan_coeffs = zeros(numahu,numcoeffs)
	chiller_coeffs = zeros(numahu)
	# Pm = zeros(numzones)
	for k = 1:numcoeffs
		for z = 1:numzones
			coeffs[z,k] = CSV.read(Coeffs_file1, DataFrame)[!,Symbol("$z")][k];
		end
	end
	for f = 1:numahu
		for k = 1:numcoeffs
			fan_coeffs[f,k] = CSV.read(Coeffs_file2, DataFrame)[!,Symbol("AHU$f")][k];
		end
		chiller_coeffs[f] = CSV.read(Coeffs_file3, DataFrame)[!,Symbol("AHU$f")][1];
	end
	# # Tinit = repeat(Tinit, Int(numzones/10));
	# global u_t = zeros(numzones)
	# global Opcost = zeros(1,Ntime)
	global Mf_p = zeros(numzones)
	global Tset_p = zeros(numzones)
	# global cputime = 0;

	global Optimization_input = Dict("numzones" => numzones,
									"numahu"    =>numahu,
									"ahuzonelist"   =>ahuzonelist,
									"coeffs"    =>coeffs,
									"fan_params"=>fan_coeffs,
									"chiller_params"=>chiller_coeffs,
									"price" 	=> price,
									"min_idx"	=>min_idx,
									# "hour_idx"=>hour_idx,
									"dT" 	=> dT,
									"t_0"	=> t,
									"H"		=> H,
									"T_init"=>Tinit,
									"T_oa"  => T_oa,
									"T_oa_p"=> T_oa_p,
									"Tzmin" => Tzmin,
									"Tzmax" => Tzmax,
									"dischargetemp" =>Tdischarge,
									"damper_min" 	=> 0,
									"damper_max" 	=> 1,
									"zoneflow_min" 	=> 0,
									"zoneflow_max" 	=> 1.3,
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
