## Objective cost function
# expression for sum of mass flows at each stage
@NLexpression(m, sum_zoneflows[f = 1:p.numfloors, t = 1:o.numstages], sum(zoneflow[f, z, window_index[t]] for z in 1:p.numzones))

# expression for fan energy (= fan power x sampling interval) consumption
@NLexpression(m, fanenergy[f = 1:p.numfloors, t = 1:o.numstages],
                (p.fan_params[f,1] + p.fan_params[f,2] * sum_zoneflows[f, t] + p.fan_params[f,3] * (sum_zoneflows[f, t])^2 + p.fan_params[f,4] * ahupressures[f, t]))

# expression for chiller capacity (energy) delivered
@NLexpression(m, chillercapacity[f = 1:p.numfloors, t = 1:o.numstages],
                p.specheat * sum_zoneflows[f, t] * (ahumixedtemp[f, t] - ahusupplytemp[f, window_index[t]]))
# # expression for  final nonnegative chiller capacity (energy) delivered
# @variable(m, chillercap_noneg[f = 1:p.numfloors, t = 1:o.numstages]>= 0.0)
# @NLconstraint(m, chillercap_nonneg_con[f = 1:p.numfloors, t = 1:o.numstages], chillercap_noneg[f,t] - chillercapacity[f,t] <= 0.0)

# expression for VAV-box reheat (energy) consumption
@NLexpression(m, vavcapacity[f = 1:p.numfloors, t = 1:o.numstages],
                p.specheat * sum(zoneflow[f, z, window_index[t]] * (zonedischargetemp[f, z, window_index[t]] - ahusupplytemp[f, window_index[t]]) for z in 1:p.numzones))

# expression for total heating capacity (AHU + VAV) (energy) delivered
@NLexpression(m, heatingcapacity[f = 1:p.numfloors, t = 1:o.numstages],
                vavcapacity[f, t])
# # chiller capacity is nonnegative
# @constraint(m, chillercap_nonneg[f = 1:p.numfloors, t = 1:o.numstages], chillercap_uppbound[f,t] - max_nonneg(chillercapacity[f,t]) <= 0)

# # expression for total chiller capacity per stage
# @NLexpression(m, chillercap_total[t = 1:o.numstages],
#                 sum(chillercapacity[f, t] for f = 1:p.numfloors))
# # expression for total chiller power per stage
# @NLexpression(m, chillerpower[t = 1:o.numstages],
#                 p.intercept + p.linear * chillercap_total[t] + p.quad * chillercap_total[t] ^ 2
#                 + p.cubic * chillercap_total[t] ^ 3)
# max_nonneg(x) = max(x, 0)
# JuMP.register(m, :max_nonneg, 1, max_nonneg, autodiff=true)

# set model objetive for different type of models
# Model 1 minimizes HVAC energy capacity + fan energy +  comfort deviation
# @NLobjective(m, Min, 60.0 * price * sum(o.weight_heat * heatingcapacity[f, t] + o.weight_cool * chillercap_noneg[f, t]+ o.weight_fan  * fanenergy[f, t] for f = 1:p.numfloors, t = 1:o.numstages) + o.penalty * slack^2)
# @NLobjective(m, Min, 60.0 * price * sum(o.weight_heat * heatingcapacity[f, t] + o.weight_cool * max_nonneg(chillercapacity[f, t]) + o.weight_fan  * fanenergy[f, t] for f = 1:p.numfloors, t = 1:o.numstages) + o.penalty * slack^2)
# @NLobjective(m, Min, 60.0 * price * sum(o.weight_heat * heatingcapacity[f, t] + o.weight_cool * chillercapacity[f, t]+ o.weight_fan * fanenergy[f, t] for f = 1:p.numfloors, t = 1:o.numstages) + o.penalty * slack^2)
@NLobjective(m, Min, 60.0 * price * sum(o.weight_heat * heatingcapacity[f, t] + o.weight_cool * (38.3671001976627+0.0150602075109015*chillercapacity[f, t]+0.000126646740696293*chillercapacity[f, t]^2 -2.72597657245502e-08*chillercapacity[f, t]^3)+ o.weight_fan * fanenergy[f, t] for f = 1:p.numfloors, t = 1:o.numstages) +
                o.penalty * sum(slack[f,z,t]^2 for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages) + 1*sum(zonedischargetemp[f, z, window_index[t]] - ahusupplytemp[f, window_index[t]] for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages))
