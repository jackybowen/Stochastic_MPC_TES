# using JuMP
# using Ipopt
# using Printf

function ModelVB!(m::Model, o::LoadParams.MPCParams, p::LoadParams.BldgParams)
## list of supervisory setpoints set by MPC (these are sent to Modelica)
    # supply-air temperatures at the AHU level
    @variable(m, ahusupplytemp[f = 1:p.numfloors, h = 1:o.numwindows],
                lower_bound = p.ahusupplytemp_min, upper_bound = p.ahusupplytemp_max)

    # damper positions at the AHU level
    @variable(m, ahudamper[f = 1:p.numfloors, h = 1:o.numwindows],
                lower_bound = p.damper_min[f], upper_bound = p.damper_max[f])

    # discharge-air temperatures at the zone level
    @variable(m, zonedischargetemp[f = 1:p.numfloors, z = 1:p.numzones, h = 1:o.numwindows],
                 upper_bound = p.zonedischargetemp_max)

    # mass-flow rates at the zone level
    @variable(m, zoneflow[f = 1:p.numfloors, z = 1:p.numzones, h = 1:o.numwindows],
                 lower_bound = p.zoneflow_min[f, z], upper_bound = p.zoneflow_max[f, z])

    # # static pressures at the AHU level
    # # note static pressure is not a true decision variable; it is computed using convex hull info once MPC is solved
    # @NLparameter(m, ahupressures[f = 1:p.numfloors, t = 1:o.numstages] == 0.0)

## list of  auxiliary optimization variables (these are NOT sent to Modelica)
    # mixed-air temperatures at the AHU level
    @variable(m, ahumixedtemp[f = 1:p.numfloors, t = 1:o.numstages])

    # return-air temperatures at the AHU level
    @variable(m, ahureturntemp[f = 1:p.numfloors, t = 1:o.numstages])

    # # internal-air temperatures at the AHU level
    # @variable(m, ahuinternaltemp[f = 1:p.numfloors, t = 1:o.numstages])

    # # Upper bound of chiller power capacity
    # @variable(m, chillercap_uppbound[f = 1:p.numfloors, t = 1:o.numstages])
    # zone temperatures
    @variable(m, zonetemp[f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages])

    # slack variable
    @variable(m, slack[f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages] >= 0.0)


    ## initialize model parameters with garbage values
    # Energy price
    @NLparameter(m, price == 1.0)

    # # Hour_of_the_day index
    # @NLparameter(m, hour_index[t = 1:o.numstages] == 1)

    # outside-air temperature
    @NLparameter(m, oat[t = 1:o.numstages] == 0.0)

    # internal loads
    @NLparameter(m, intload[f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages] == 0.0)

    # initial zone temperature
    @NLparameter(m, zonetemp_init[f = 1:p.numfloors, z = 1:p.numzones] == 0.0)
    # initial zone flow
    @NLparameter(m, zoneflow_init[f = 1:p.numfloors, z = 1:p.numzones] == 0.0)
    # initial zone discharge air temperature
    @NLparameter(m, zonedischargetemp_init[f = 1:p.numfloors, z = 1:p.numzones] == 0.0)

    # default min and max zone temperatures
    @NLparameter(m, heatsetpoint[t = 1:o.numstages] == 0.0)
    @NLparameter(m, coolsetpoint[t = 1:o.numstages] == 0.0)
## Constraints
    # zone temperature constraints in stage 1
    @NLconstraint(m, zone_cons_first[f = 1:p.numfloors, z = 1:p.numzones, t = 1:1],
                    expr[f, z, t] == p.zonetemp_alpha[f, z] * zonetemp_init[f, z])

    # zone temperature constraints in other stages
    @NLconstraint(m, zone_cons[f = 1:p.numfloors, z = 1:p.numzones, t = 2:o.numstages],
                    expr[f, z, t] == p.zonetemp_alpha[f, z] * zonetemp[f, z, t-1])

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

## Objective term
    # expression for sum of mass flows at each stage
    @NLexpression(m, sum_zoneflows[f = 1:p.numfloors, t = 1:o.numstages], sum(zoneflow[f, z, window_index[t]] for z in 1:p.numzones))

    # expression for fan energy (= fan power x sampling interval) consumption
    @NLexpression(m, fanenergy[f = 1:p.numfloors, t = 1:o.numstages],
                    (p.fan_params[f,1] + p.fan_params[f,2] * sum_zoneflows[f, t] + p.fan_params[f,3] * (sum_zoneflows[f, t])^2 + p.fan_params[f,4] * (sum_zoneflows[f, t])^3))

    # expression for chiller capacity (energy) delivered
    @NLexpression(m, chillercapacity[f = 1:p.numfloors, t = 1:o.numstages],
                    p.specheat * sum_zoneflows[f, t] * (ahumixedtemp[f, t] - ahusupplytemp[f, window_index[t]])/(eta_cop))
    # # expression for  final nonnegative chiller capacity (energy) delivered
    # @variable(m, chillercap_noneg[f = 1:p.numfloors, t = 1:o.numstages]>= 0.0)
    # @NLconstraint(m, chillercap_nonneg_con[f = 1:p.numfloors, t = 1:o.numstages], chillercap_noneg[f,t] - chillercapacity[f,t] <= 0.0)
    #
    # # expression for VAV-box reheat (energy) consumption
    # @NLexpression(m, vavcapacity[f = 1:p.numfloors, t = 1:o.numstages],
    #                 p.specheat * sum(zoneflow[f, z, window_index[t]] * (zonedischargetemp[f, z, window_index[t]] - ahusupplytemp[f, window_index[t]]) for z in 1:p.numzones))
    #
    # # expression for total heating capacity (AHU + VAV) (energy) delivered
    # @NLexpression(m, heatingcapacity[f = 1:p.numfloors, t = 1:o.numstages],vavcapacity[f, t])
    # @NLexpression(m, chillercapacity,(38.3671001976627+0.0150602075109015*chillercapacity[f, t]+0.000126646740696293*chillercapacity[f, t]^2 -2.72597657245502e-08*chillercapacity[f, t]^3)+ o.weight_fan * fanenergy[f, t] for f = 1:p.numfloors, t = 1:o.numstages))
    return m
end
