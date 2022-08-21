module Opt_engine

# using Clp
using Ipopt
# using ECOS
using JuMP
using DataFrames
using CSV
# using Conda
# using PyCall
using Statistics
# Conda.add("nomkl")
# Conda.add("joblib")
# Conda.add("scikit-learn")
# Conda.rm("mkl")
cd(dirname(@__FILE__))



function run(GUI_input::Dict)

    ###### read inputs ######

    for (key, value) in GUI_input
        @eval $(Symbol(key)) = $value
    end

    ###### Chiller power Optimal Scheduling ######
    # looking at a 24-hour-ahead time window for the optimal scheduling
    t0 = Int(t_0); # start timestamp in minute
    if solverflag == 1
        m = Model(Ipopt.Optimizer);
        set_optimizer_attribute(m, "print_level", 3)
        set_optimizer_attribute(m, "max_iter", 1000)
		# set_optimizer_attribute(m, "linear_solver", "ma86")
        # set_optimizer_attribute(m, "threads", 4)
    	set_time_limit_sec(m, 300.0)
    elseif solverflag == 2
        m = Model(Clp.Optimizer);
        set_optimizer_attribute(m, "LogLevel", 1);
    end
    # m = Model(Clp.Optimizer);
    # set_optimizer_attribute(m, "LogLevel", 0);

	# T_noise = zeros(H)
	## list of supervisory setpoints set by MPC (these are sent to Modelica)
	# supply-air temperatures at the AHU level
	# @variable(m, ahusupplytemp[t = 1:H],lower_bound = ahusupplytemp_min, upper_bound = ahusupplytemp_max)

	# damper positions at the AHU level
	@variable(m, ahudamper[f = 1:numahu, t = 1:H], lower_bound = damper_min, upper_bound = damper_max)

	# discharge-air temperatures at the zone level
	@variable(m, ahudischargetemp[f = 1:numahu, t = 1:H])

	# mass-flow rates at the zone level
	@variable(m, zoneflow[z = 1:numzones, t = 1:H],lower_bound = zoneflow_min, upper_bound = zoneflow_max)

	## list of  auxiliary optimization variables (these are NOT sent to Modelica)
	# mixed-air temperatures at the AHU level
	@variable(m, ahumixedtemp[f = 1:numahu, t = 1:H])

	# return-air temperatures at the AHU level
	@variable(m, ahureturntemp[f = 1:numahu, t = 1:H])

    # @variable(m, p_c[z = 1:numzones, t = 1:H], Bin);   # hourly chiller power (Binary variable)
	@variable(m, slack_u[z = 1:numzones, t = 1:H] >= 0.0)
    @variable(m, slack_d[z = 1:numzones, t = 1:H] >= 0.0)
    @variable(m, Tzon[z = 1:numzones,t = 1:H])#,lower_bound = 16,upper_bound = 26)   # hourly setpoint temperatures(unit in C)
	@constraint(m, maxcomfort_cons[z=1:numzones, t=1:H], Tzon[z,t] <= Tzmax[z,t]+slack_u[z,t])
    @constraint(m, mincomfort_cons[z=1:numzones, t=1:H], Tzon[z,t] >= Tzmin[z,t]-slack_d[z,t])

	# if  uncertout_flag == 0
    # @constraint(m, zon_temp_cons_init[z=1:numzones, t=1], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*T_init[z]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t])
    # @constraint(m, zon_temp_cons[z=1:numzones, t=2:H], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*Tzon[z,t-1]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t])

	@constraint(m, mflow_cons_init[z in mzonelist, t=1], zoneflow[z,t] == coeffs[z,4] + coeffs[z,1]*T_init[z]+coeffs[z,2]*T_oa_p[t]+ coeffs[z,3]*Tzon[z,t])
	@constraint(m, mflow_cons[z in mzonelist, t=2:H], zoneflow[z,t] == coeffs[z,4] + coeffs[z,1]*Tzon[z,t-1]+coeffs[z,2]*T_oa_p[t] + coeffs[z,3]*Tzon[z,t])
	# Constant air mass flow rate for single zone(AHU002, AHU004)
	@constraint(m, mflow_cons_const[z in szonelist, t=1:H], zoneflow[z,t] == m_init[z])

	# Fan power constraints
	@NLexpression(m, sum_zoneflows[f=1:numahu, t = 1:H], sum(zoneflow[z, t] for z in ahuzonelist[f]))
	@NLexpression(m, fanenergy[f=1:numahu, t = 1:H], fan_params[f,4] + fan_params[f,1] * sum_zoneflows[f,t] + fan_params[f,2] * (sum_zoneflows[f,t])^2 + fan_params[f,3] * (sum_zoneflows[f,t])^3)

	# Chiller power constraints
    # return-air temperature constraints (definiton)
    @NLconstraint(m, return_cons[f = 1:numahu, t = 1:H],ahureturntemp[f,t] * sum_zoneflows[f,t] == sum(zoneflow[z,t] * Tzon[z,t] for z in ahuzonelist[f]))
    # mixed-air temperature constraints (definition)
    @NLconstraint(m, mixed_cons[f = 1:numahu, t = 1:H],ahumixedtemp[f,t] == ahudamper[f,t] * T_oa_p[t] + (1.0 - ahudamper[f,t]) * ahureturntemp[f,t])
	# discharge-air temperature constraints for single-zone served AHU(AHU-002, AHU-004)
	@constraint(m, supplytemp_cons_init[f in sahulist, t=1], ahudischargetemp[f,t] == coeffs[ahuzonelist[f],4] + coeffs[ahuzonelist[f],1]*T_init[ahuzonelist[f]]+coeffs[ahuzonelist[f],2]*T_oa_p[t] + coeffs[ahuzonelist[f],3]*Tzon[ahuzonelist[f],t])
	@constraint(m, supplytemp_cons[f in sahulist, t=2:H], ahudischargetemp[f,t] == coeffs[ahuzonelist[f],4] + coeffs[ahuzonelist[f],1]*Tzon[ahuzonelist[f],t-1]+coeffs[ahuzonelist[f],2]*T_oa_p[t] + coeffs[ahuzonelist[f],3]*Tzon[ahuzonelist[f],t])
	# Baseline discharge-air temperature for multi-zone served AHU(AHU001, AHU003)
	@constraint(m, supplytemp_cons_other[f in mahulist, t=1:H], ahudischargetemp[f,t] == dischargetemp[f,t])

	@constraint(m, positive_load_bound[f=1:numahu,t=1:H], ahumixedtemp[f,t] - ahudischargetemp[f,t] >= 0) # if Tmixed - Tdis < 0 , let load = 0

	# expression for AHU load demand
	@NLexpression(m, load[t = 1:H], sum(sum_zoneflows[f,t] * (ahumixedtemp[f,t] - ahudischargetemp[f,t]) for f in 1:numahu))
	# expression for chiller capacity (energy) delivered
	@NLexpression(m, chillercapacity[t = 1:H], chiller_params[1] * load[t] + chiller_params[2] * load[t]^2 + chiller_params[3])

	# elseif uncertout_flag == 1
	# 	@constraint(m, zon_temp_cons_init[z=1:numzones, t=1], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*T_init[z]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t] + Tz_n[z,t])
    #     @constraint(m, zon_temp_cons[z=1:numzones, t=2:H], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*Tzon[z,t-1]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t] + Tz_n[z,t])
	# end

    # # Force heat power to be zero
    # @constraint(m, heater_cons[z = 1:numzones, t=1:H], p_h[z,t] == 0.0)
	penalty_u = 1e8;
    penalty_d = 1e7;
    if optcost_flag == 1
		# for t=1:H
		# 	println("Time:",t0+Int((t-1)*60*dT))
		# end
        @expression(m, pr[t=1:H], price[min_idx[t0+Int((t-1)*60*dT)+1]])
        @NLobjective(m, Min, sum(pr[t]*(fanenergy[f,t] + chillercapacity[t]) for f=1:numahu, t=1:H) + sum(penalty_u*slack_u[z,t]+penalty_d*slack_d[z,t] for z = 1:numzones, t = 1:H));
        # @objective(m, Min, sum((p_c[z,t]+p_h[z,t]) for z=1:numzones, t=1:H));
    # else
    #     @objective(m, Min, sum((Pm[z]*p_c[z,t]) for z=1:numzones, t=1:H)+penalty*sum(slack[z,t]^2 for z = 1:numzones, t = 1:H));
    end

    status = optimize!(m);
    println("Optimal scheduling for minute ", Int(t0), " -- Status: ", termination_status(m));
    sol_mflow = JuMP.value.(zoneflow); # println(sol_mflow);
    sol_tzon = JuMP.value.(Tzon); # println(sol_tzon);
    sol_pr = JuMP.value.(pr);
    sol_slack_u = JuMP.value.(slack_u);
    sol_slack_d = JuMP.value.(slack_d);
	sol_Tz_upp = Tzmax.+sol_slack_u
	sol_Tz_low = Tzmin.-sol_slack_d
	# println("Upper bound:",sol_Tz_upp)
	# println("Lower bound", sol_Tz_low)

	# Nmin = Int(dT_m/dT);
	# sol_mflow = zeros(numzones,Hm);
	# sol_tzon = zeros(numzones,Hm);
	# sol_pr = mean(reshape(sol_pr_all,Nmin,:),dims=1);
	# sol_slack_u = zeros(numzones,Hm);
	# sol_slack_d = zeros(numzones,Hm);
	# sol_Tz_upp = zeros(numzones,Hm);
	# for z = 1:numzones
	# 	sol_mflow[z,:] = mean(reshape(sol_mflow_all[z,:],Nmin,:),dims=1); # println(sol_mflow);
	#     sol_tzon[z,:] = mean(reshape(sol_tzon_all[z,:],Nmin,:),dims=1);# println(sol_tzon);
	#     sol_slack_u[z,:] = mean(reshape(sol_slack_u_all[z,:],Nmin,:),dims=1);
	#     sol_slack_d[z,:] = mean(reshape(sol_slack_d_all[z,:],Nmin,:),dims=1);
	# 	sol_Tz_upp[z,:] = mean(reshape(sol_Tz_upp_all[z,:],Nmin,:),dims=1);
	# end
    sol_fanpower = JuMP.value.(fanenergy)
	sol_sumflow = JuMP.value.(sum_zoneflows)
	# println("Fanpower cost 01: ", sum(sol_fanpower[1,t] for t=1:H))
	# println("Fanpower cost 02: ", sum(sol_fanpower[2,t] for t=1:H))
	# println("Fanpower cost 03: ", sum(sol_fanpower[3,t] for t=1:H))
	# println("Fanpower cost 04: ", sum(sol_fanpower[4,t] for t=1:H))
	println("Price: ", sol_pr)
	println("SumZoneFlow 01: ", sol_sumflow[1,:])
	println("SumZoneFlow 02: ", sol_sumflow[2,:])
	println("SumZoneFlow 03: ", sol_sumflow[3,:])
	println("SumZoneFlow 04: ", sol_sumflow[4,:])
	# println("Fan 003 parameter intercept: ",fan_params[3,4])
	# println("Fan 003 parameter 1st order: ",fan_params[3,1] * sol_sumflow[3,:])
	# println("Fan 003 parameter 2nd order: ",fan_params[3,2] * (sol_sumflow[3,:]).^2)
	# println("Fan 003 parameter 3rd order: ",fan_params[3,3] * (sol_sumflow[3,:]).^3)

	sol_dischargetemp = JuMP.value.(ahudischargetemp)
	sol_mixedtemp = JuMP.value.(ahumixedtemp)
	println("AHU 001 discharge temperature: ",sol_dischargetemp[1,:])
	println("AHU 002 discharge temperature: ",sol_dischargetemp[2,:])
	println("AHU 003 discharge temperature: ",sol_dischargetemp[3,:])
	println("AHU 004 discharge temperature: ",sol_dischargetemp[4,:])
	println("AHU 001 mixed temperature: ",sol_mixedtemp[1,:])
	println("AHU 002 mixed temperature: ",sol_mixedtemp[2,:])
	println("AHU 003 mixed temperature: ",sol_mixedtemp[3,:])
	println("AHU 004 mixed temperature: ",sol_mixedtemp[4,:])
	sol_load = JuMP.value.(load)
	println("Load: ", sol_load)
	sol_chillerpower=JuMP.value.(chillercapacity)
	println("Chillerpower cost:", sol_chillerpower)
    sol_cost = sum(sol_pr[t]*(sol_fanpower[f,t]+sol_chillerpower[t]) for f=1:numahu, t=1:H)
    sol_cost1 = sum(penalty_u*sol_slack_u[z,t] for z=1:numzones,t=1:H)
    sol_cost2 = sum(penalty_d*sol_slack_d[z,t] for z=1:numzones,t=1:H)
    println("Total Power cost: ",sol_cost)
    println("Upper bound relaxation cost: ",sol_cost1)
    println("Lower bound relaxation cost: ",sol_cost2)

    df_cost = DataFrame(t0 = t0, cost = sol_cost, cost1 = sol_cost1, cost2 = sol_cost2)
    CSV.write("saved_cost.csv", df_cost, append=true)
    # result = DataFrame(

    result1 = DataFrame()
    result2 = DataFrame()
    result3 = DataFrame()
    result4 = DataFrame()
    result5 = DataFrame()
    result6 = DataFrame()
    result7 = DataFrame()
	result8 = DataFrame()
    default_Values = 0;
    for z = 1:Int(numzones)
        insertcols!(result1,z,Symbol("sol_mflow_$z") =>vcat(m_init[z],sol_mflow[z,:]))
        insertcols!(result2,z,Symbol("sol_tzon_$z")=>vcat(T_init[z],sol_tzon[z,:]))
        insertcols!(result3,z,Symbol("sol_slack_u_$z")=>vcat(default_Values,sol_slack_u[z,:]))
        insertcols!(result4,z,Symbol("sol_slack_d_$z")=>vcat(default_Values,sol_slack_d[z,:]))
		insertcols!(result5,z,Symbol("upp_penalty_cost_$z")=>vcat(default_Values,penalty_u*sol_slack_u[z,:]))
        insertcols!(result6,z,Symbol("low_penalty_cost_$z")=>vcat(default_Values,penalty_d*sol_slack_d[z,:]))
        insertcols!(result7,z,Symbol("Tupp_$z")=>vcat(Tzmax[z,1],sol_Tz_upp[z,:]))
        insertcols!(result8,z,Symbol("Tlow_$z")=>vcat(Tzmin[z,1],sol_Tz_low[z,:]))
    end
    result0 = DataFrame(price = vcat(price[min_idx[t0]], vec(sol_pr)), OutdoorTemperature = vcat(T_oa,vec(T_oa_p)))
    result = hcat(result0,result1,result2,result3,result4,result5,result6,result7,result8)
	# println(T_noise)
    return result;

end
end
