module Opt_engine

using Clp
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
    fan_params = [0 0.9468 0.0186 0]
    # looking at a 24-hour-ahead time window for the optimal scheduling
    H = Int(24/dT);
    t0 = Int(t_0); # start timestamp in minute
    if solverflag == 1
        m = Model(Ipopt.Optimizer);
        set_optimizer_attribute(m, "print_level", 3)
        set_optimizer_attribute(m, "max_iter", 1000)
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
	@variable(m, ahudamper[t = 1:H], lower_bound = damper_min, upper_bound = damper_max)

	# # discharge-air temperatures at the zone level
	# @variable(m, zonedischargetemp[z = 1:numzones, t = 1:H], upper_bound = zonedischargetemp_max)

	# mass-flow rates at the zone level
	@variable(m, zoneflow[z = 1:numzones, t = 1:H],lower_bound = zoneflow_min, upper_bound = zoneflow_max)

	## list of  auxiliary optimization variables (these are NOT sent to Modelica)
	# mixed-air temperatures at the AHU level
	@variable(m, ahumixedtemp[t = 1:H])

	# return-air temperatures at the AHU level
	@variable(m, ahureturntemp[t = 1:H])

    # @variable(m, p_c[z = 1:numzones, t = 1:H], Bin);   # hourly chiller power (Binary variable)
	# @variable(m, slack_u[z = 1:numzones, t = 1:H] >= 0.0)
    # @variable(m, slack_d[z = 1:numzones, t = 1:H] >= 0.0)
    @variable(m, Tzon[z = 1:numzones,t = 1:H])#,lower_bound = 16,upper_bound = 26)   # hourly setpoint temperatures(unit in C)
	@constraint(m, maxcomfort_cons[z=1:numzones, t=1:H], Tzon[z,t] <= Tzmax[z,t])
    @constraint(m, mincomfort_cons[z=1:numzones, t=1:H], Tzon[z,t] >= Tzmin[z,t])

	# if  uncertout_flag == 0
    # @constraint(m, zon_temp_cons_init[z=1:numzones, t=1], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*T_init[z]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t])
    # @constraint(m, zon_temp_cons[z=1:numzones, t=2:H], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*Tzon[z,t-1]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t])
	@constraint(m, mflow_cons_init[z=1:numzones, t=1], zoneflow[z,t] == coeffs[z,4] + coeffs[z,1]*T_init[z]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,3]*Tzon[z,t])
	@constraint(m, mflow_cons[z=1:numzones, t=2:H], zoneflow[z,t] == coeffs[z,4] + coeffs[z,1]*Tzon[z,t-1]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,3]*Tzon[z,t])
	# VB constraints
	@NLexpression(m, sum_zoneflows[t = 1:H], sum(zoneflow[z, t] for z in 1:numzones))
	@NLexpression(m, fanenergy[t = 1:H], fan_params[1] + fan_params[2] * sum_zoneflows[t] + fan_params[3] * (sum_zoneflows[t])^2 + fan_params[4] * (sum_zoneflows[t])^3)

    # return-air temperature constraints (definiton)
    @NLconstraint(m, return_cons[t = 1:H],ahureturntemp[t] * sum(zoneflow[z,t] for z in 1:numzones) == sum(zoneflow[z,t] * Tzon[z,t] for z in 1:numzones))
    # mixed-air temperature constraints (definition)
    @NLconstraint(m, mixed_cons[t = 1:H],ahumixedtemp[t] == ahudamper[t] * T_oa_p[t] + (1.0 - ahudamper[t]) * ahureturntemp[t])
	# expression for chiller capacity (energy) delivered
	@NLexpression(m, chillercapacity[t = 1:H], specheat/COP * sum_zoneflows[t] * (ahumixedtemp[t] - 12.8))

	# elseif uncertout_flag == 1
	# 	@constraint(m, zon_temp_cons_init[z=1:numzones, t=1], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*T_init[z]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t] + Tz_n[z,t])
    #     @constraint(m, zon_temp_cons[z=1:numzones, t=2:H], Tzon[z,t] == coeffs[z,4] + coeffs[z,3]*Tzon[z,t-1]+coeffs[z,2]*(T_oa_p[t]) + coeffs[z,1]*Pm[z]*p_c[z,t] + Tz_n[z,t])
	# end

    # # Force heat power to be zero
    # @constraint(m, heater_cons[z = 1:numzones, t=1:H], p_h[z,t] == 0.0)
	# penalty_u = 1e5;
    # penalty_d = 1e3;
    if optcost_flag == 1
		# for t=1:H
		# 	println("Time:",t0+Int((t-1)*60*dT))
		# end
        @expression(m, pr[t=1:H], price[min_idx[t0+Int((t-1)*60*dT)+1]])
        @NLobjective(m, Min, sum(pr[t]*(fanenergy[t] + chillercapacity[t]) for t=1:H))#+sum(penalty_u*slack_u[z,t]+penalty_d*slack_d[z,t] for z = 1:numzones, t = 1:H));
        # @objective(m, Min, sum((p_c[z,t]+p_h[z,t]) for z=1:numzones, t=1:H));
    # else
    #     @objective(m, Min, sum((Pm[z]*p_c[z,t]) for z=1:numzones, t=1:H)+penalty*sum(slack[z,t]^2 for z = 1:numzones, t = 1:H));
    end

    status = optimize!(m);
    println("Optimal scheduling for minute ", Int(t0), " -- Status: ", termination_status(m));
    sol_mflow = JuMP.value.(zoneflow); # println(sol_mflow);
    sol_tzon = JuMP.value.(Tzon); # println(sol_tzon);
    sol_pr = JuMP.value.(pr);
    # sol_slack_u = JuMP.value.(slack_u);
    # sol_slack_d = JuMP.value.(slack_d);
	sol_Tz_upp = Tzmax
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
    # sol_cost = sum(sol_pr.*(Pm[z]*sol_mflow[z,:]) for z=1:numzones)
    # sol_cost = sum(sol_pr[t]*(Pm[z]*sol_mflow[z,t]) for z=1:numzones,t=1:Hm)
    # sol_cost1 = sum(sol_slack_u[z,t] for z=1:numzones, t=1:Hm)
    # sol_cost2 = sum(sol_slack_d[z,t] for z = 1:numzones, t = 1:Hm)
    # println(sol_cost)
    # println(sol_cost1)
    # println(sol_cost2)

	sol_Tz_low = Tzmin
    # df_cost = DataFrame(t0 = t0, cost = sol_cost, cost1 = sol_cost1, cost2 = sol_cost2)
    # CSV.write("saved_cost.csv", df_cost, append=true)
    # result = DataFrame(
    # sol_mflow_1 = sol_mflow[1,:],
    # sol_mflow_2 = sol_mflow[2,:],
    # sol_mflow_3 = sol_mflow[3,:],
    # sol_mflow_4 = sol_mflow[4,:],
    # sol_mflow_5 = sol_mflow[5,:],
    # sol_mflow_6 = sol_mflow[6,:],
    # sol_mflow_7 = sol_mflow[7,:],
    # sol_mflow_8 = sol_mflow[8,:],
    # sol_mflow_9 = sol_mflow[9,:],
    # sol_mflow_10= sol_mflow[10,:],
    # # sol_p_h_1 = POW_h[1,:],
    # # sol_p_h_2 = POW_h[2,:],
    # # sol_p_h_3 = POW_h[3,:],
    # # sol_p_h_4 = POW_h[4,:],
    # # sol_p_h_5 = POW_h[5,:],
    # # sol_p_h_6 = POW_h[6,:],
    # # sol_p_h_7 = POW_h[7,:],
    # # sol_p_h_8 = POW_h[8,:],
    # # sol_p_h_9 = POW_h[9,:],
    # # sol_p_h_10 = POW_h[10,:],
    # sol_tzon_1 = sol_tzon[1,:],
    # sol_tzon_2 = sol_tzon[2,:],
    # sol_tzon_3 = sol_tzon[3,:],
    # sol_tzon_4 = sol_tzon[4,:],
    # sol_tzon_5 = sol_tzon[5,:],
    # sol_tzon_6 = sol_tzon[6,:],
    # sol_tzon_7 = sol_tzon[7,:],
    # sol_tzon_8 = sol_tzon[8,:],
    # sol_tzon_9 = sol_tzon[9,:],
    # sol_tzon_10= sol_tzon[10,:]);
    result1 = DataFrame()
    result2 = DataFrame()
    result3 = DataFrame()
    result4 = DataFrame()
    result5 = DataFrame()
    result6 = DataFrame()
    result7 = DataFrame()
    default_Values = -100000;
    for z = 1:Int(numzones)
        insertcols!(result1,z,Symbol("sol_mflow_$z") =>vcat(default_Values,sol_mflow[z,:]))
        insertcols!(result2,z,Symbol("sol_tzon_$z")=>vcat(T_init[z],sol_tzon[z,:]))
        # insertcols!(result3,z,Symbol("sol_slack_u_$z")=>vcat(default_Values,sol_slack_u[z,:]))
        # insertcols!(result4,z,Symbol("sol_slack_d_$z")=>vcat(default_Values,sol_slack_d[z,:]))
        insertcols!(result5,z,Symbol("Tsp_$z") =>vcat(default_Values,Tzmax[z,:]))
        insertcols!(result6,z,Symbol("Tupp_$z")=>vcat(default_Values,sol_Tz_upp[z,:]))
        insertcols!(result7,z,Symbol("Tlow_$z")=>vcat(default_Values,sol_Tz_low[z,:]))
    end
    result0 = DataFrame(price = vcat(price[min_idx[t0]], vec(sol_pr)), OutdoorTemperature = vcat(T_oa,vec(T_oa_p)))
    result = hcat(result0,result1,result2,result5,result6,result7)
	# println(T_noise)
    return result;

end
end
