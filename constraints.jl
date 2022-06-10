## model constraints
include("./zone_temp_constraints.jl")
# # zone temperature dynamic constraints
# expr = Dict{Tuple{Int64, Int64, Int64}, JuMP.NonlinearExpression}()  # initialize an empty expression
# for f in 1:p.numfloors, z in 1:p.numzones, t in 1:o.numstages
#     h = window_index[t]
#     # expression for zone dynamics minus the alpha * (zonetemp[t-1]) term
#     expr[f, z, t]  = @NLexpression(m,
#                     zonetemp[f, z, t] - p.zonetemp_beta[f, z] * zoneflow[f, z, h] * (zonedischargetemp[f, z, h] - zonetemp[f, z, t])
#                     - p.zonetemp_gamma[f, z] * oat[t] - intload[f, z, t])
# end
#
# # zone temperature constraints in stage 1
# @NLconstraint(m, zone_cons_first[f = 1:p.numfloors, z = 1:p.numzones, t = 1:1],
#                 expr[f, z, t] == p.zonetemp_alpha[f, z] * zonetemp_init[f, z])
#
# # zone temperature constraints in other stages
# @NLconstraint(m, zone_cons[f = 1:p.numfloors, z = 1:p.numzones, t = 2:o.numstages],
#                 expr[f, z, t] == p.zonetemp_alpha[f, z] * zonetemp[f, z, t-1])

# return-air temperature constraints (definiton)
@NLconstraint(m, return_cons[f = 1:p.numfloors, t = 1:o.numstages],
                ahureturntemp[f, t] * sum(zoneflow[f, z, window_index[t]] for z in 1:p.numzones) == sum(zoneflow[f, z, window_index[t]] * zonetemp[f, z, t] for z in 1:p.numzones))

# mixed-air temperature constraints (definition)
@NLconstraint(m, mixed_cons[f = 1:p.numfloors, t = 1:o.numstages],
                ahumixedtemp[f, t] == ahudamper[f, window_index[t]] * oat[t] + (1.0 - ahudamper[f, window_index[t]]) * ahureturntemp[f, t])

# comfort bounding constraints
@NLconstraint(m, mincomfort_cons[f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages] ,zonetemp[f, z, t] >= heatsetpoint[t] - slack[f,z,t])
@NLconstraint(m, maxcomfort_cons[f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages], zonetemp[f, z, t] <= coolsetpoint[t] + slack[f,z,t])

# # internal temperature is upper bound of mixed-air temperature
# @constraint(m, intmixed_cons[f = 1:p.numfloors, t = 1:o.numstages], ahuinternaltemp[f, t] >= ahumixedtemp[f, t])
#
# # internal temperature is upper bound of supply-air temperature
# @constraint(m, intsupply_cons[f = 1:p.numfloors, t = 1:o.numstages], ahuinternaltemp[f, t] >= ahusupplytemp[f, window_index[t]])

# discharge temperature is upper bound of supply-air temperatures
@constraint(m, discsup_cons[f = 1:p.numfloors, z = 1:p.numzones, h = 1:o.numwindows], zonedischargetemp[f, z, h] >= ahusupplytemp[f, h] + 3)

# deviation between supply-air temp and mixed-air temp
# upper and lower bound for the difference between mixedair temp and supply-air temp
# these parameters are used to control HVAC start/stop
@NLparameter(m, bound[t = 1:o.numstages] == 0.0)
@NLconstraint(m, supmix_upper_cons[f = 1:p.numfloors, t=1:o.numstages], ahusupplytemp[f, window_index[t]] - ahumixedtemp[f, t] <= bound[t])
@NLconstraint(m, supmix_lower_cons[f = 1:p.numfloors, t=1:o.numstages], ahumixedtemp[f, t] - ahusupplytemp[f, window_index[t]] <= bound[t])

# rate of change of supply-air temperature setpoints are bounded for baseline 1 in closed loop
#if o.cl_baseline == 1
    # previously implemented supply-air temperature setpoints
    #@NLparameter(m, current_ahusupplytemp[f = 1:p.numfloors] == 0.0)
    #@NLparameter(m, true_ahusupplytemp[f = 1:p.numfloors] == max(p.ahusupplytemp_min, getvalue(current_ahusupplytemp[f])))
    ##bounds for first-stage
    #@NLconstraint(m, deltasupplytemp_lower_cons[f = 1:p.numfloors],
    #                ahusupplytemp[f, 1] - true_ahusupplytemp[f] >= -1.0 * o.cl_rate_supplytemp)
    #@NLconstraint(m, deltasupplytemp_upper_cons[f = 1:p.numfloors],
    #                ahusupplytemp[f, 1] - true_ahusupplytemp[f] <= 1.0 * o.cl_rate_supplytemp)
#end

# constraining the change between 2 consecutive AHU supply temperature set points to avoid big jumps as consequences of
# calculating for comfort or energy savings
@NLconstraint(m, [f = 1:p.numfloors, h = 1:o.numwindows-1],
            ahusupplytemp[f, h + 1] - ahusupplytemp[f, h] <= p.ahusupplytemp_max_dev)
