module Controller
using Dates
using Printf
using MathOptInterface
const MOI = MathOptInterface
# import the Julia script with the optimization model and associated dependencies
include("EngineOpt.jl")
global status_optimalSolution = [MOI.OPTIMAL]
global status_goodSolution = [MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
global status_badSolution = [MOI.ALMOST_LOCALLY_SOLVED, MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE, MOI.LOCALLY_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED, MOI.ALMOST_INFEASIBLE, MOI.ALMOST_DUAL_INFEASIBLE]
# global status_goodSolution = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL,MOI.ALMOST_LOCALLY_SOLVED]
# global status_badSolution = [MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE, MOI.LOCALLY_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED, MOI.ALMOST_INFEASIBLE, MOI.ALMOST_DUAL_INFEASIBLE]
global status_userLim = [MOI.ITERATION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.SOLUTION_LIMIT, MOI.MEMORY_LIMIT, MOI.OBJECTIVE_LIMIT, MOI.NORM_LIMIT, MOI.OTHER_LIMIT]

global originalMaxIter = o.maxiter
global cwindow = o.controlwindow
global numfloors = p.numfloors
global numzones = p.numzones
global start_minute = 0


"""
This module implements an MPC controller.
"""

function compute_control!(u::Dict, currentMeasurements::Dict)
    # compute the control input from the measurement.
    # y contains the current values of the measurements.
    # {<measurement_name>:<measurement_value>}
    # compute the control input from the measurement.
    # u: dict Defines the control input to be used for the next step.
    # {<input_name> : <input_value>}
    # global start_minute = 288001 # 17280060/60
    current_minute = currentMeasurements["Time"]/60.0

    minute = current_minute - start_minute
    minute_of_day = Helper.minute_of_day(current_minute)
    global dfCurrentMeasurements = dict2df!(dfCurrentMeasurements, currentMeasurements)
    # note df_history only requires "updated" measurements, not setpoints
    global dfCurrentMeasurements_history = minute == 0.0 ? dfCurrentMeasurements : updatehistory!(dfCurrentMeasurements_history, dfCurrentMeasurements)

    h, c =  comfortbounds(current_minute)
    pred_oat = predictambient(dfCurrentMeasurements[1, :TOutDryBul_y], minute)
    # obtain internal loads for each zone (3D array)
    pred_load = predictloads(dfCurrentMeasurements_history, minute_of_day)
    # Renew the energy price
    price = currentMeasurements["EnergyPrice"]
    # update MPC model parameters
    mpc_params = Dict("EnergyPrice"=>price, "heat_sp" => h, "cool_sp" => c, "pred_oat" => pred_oat, "pred_loads" => pred_load)
    updatemodelparams(dfCurrentMeasurements, mpc_params)
    if p.zonetempmodelflag == 1
        updatemodelinit(dfCurrentMeasurements_history, minute_of_day)
    end

    # solve MPC model
    global solverinfo = solvemodel()
    @printf("================= Quick statistics of optimization algorithms. =============\n")
    @printf("Going for %d max iterations.\n", o.maxiter)
    # println(JuMP.SimplexIterations(m))
    #println(MOI.get(m, MOI.SimplexIterations()))
    @printf("Optimization took %.4f seconds, and ended with status %s.\n", solverinfo["soltime"], solverinfo["status"])
    #@printf("size of ahusupplytemp: %d, %d, %d\n", size(ahusupplytemp, 1), size(ahusupplytemp, 2), size(ahusupplytemp, 3))
    #@printf("size of ahudamper: %d, %d, %d\n", size(ahudamper, 1), size(ahudamper, 2), size(ahudamper, 3))
    #@printf("size of zonedischargetemp: %d, %d, %d\n", size(zonedischargetemp, 1), size(zonedischargetemp, 2), size(zonedischargetemp, 3))
    #@printf("size of zoneflow: %d, %d, %d\n", size(zoneflow, 1), size(zoneflow, 2), size(zoneflow, 3))
    #@printf("size of ahupressures: %d, %d, %d\n", size(ahupressures, 1), size(ahupressures, 2), size(ahupressures, 3))
    #@printf("size of zonetemp: %d, %d, %d\n", size(zonetemp, 1), size(zonetemp, 2), size(zonetemp, 3))
    for f = 1 : p.numfloors
        @printf("floor %d -> ahusupplytemp = %.4f, constrain: [%.2f, %.2f]\n", f, JuMP.value(ahusupplytemp[f, 1]), p.ahusupplytemp_min, p.ahusupplytemp_max)
        @printf("floor %d -> ahumixedtemp = %.4f, constrain: [%.2f, %.2f]\n", f, JuMP.value(ahumixedtemp[f, 1]), JuMP.value(ahureturntemp[f, 1]), JuMP.value(oat[1]))
        @printf("floor %d -> ahudamper = %.4f, constrain: [%.2f, %2f]\n", f, JuMP.value(ahudamper[f, 1]), p.damper_min[f], p.damper_max[f])
        for z = 1 : p.numzones
            zid = (f-1)*p.numzones + z;
            if  o.debugmodeflag == 1
                df_Tinit[!,Symbol("Tinit")] = [JuMP.value(zonetemp_init[f,z])];
                CSV.write("./debuginfo/Zone$(zid)_time$(minute+1)_initT.csv",df_Tinit, append=false)

            end
            @printf("floor %d, zone %d -> zonedischargetemp = %.4f, constrain: [%.2f, %2f]\n", f, z, JuMP.value(zonedischargetemp[f, z, 1]), JuMP.value(ahusupplytemp[f, 1]), p.zonedischargetemp_max)
            @printf("floor %d, zone %d -> zoneflow = %.4f, constrain: [%.2f, %2f]\n", f, z, JuMP.value(zoneflow[f, z, 1]), p.zoneflow_min[f,z], p.zoneflow_max[f,z])
            # @printf("floor %d, zone %d -> zonetemp = %.4f,zonetemp_Real = %.4f, constrain: [%.2f, %2f]\n", f, z, JuMP.value(zonetemp[f, z, 1]), dfCurrentMeasurements[1, Symbol("floor$(f)_zon$(z)_TRooAir_y")], p.zonetemp_min_occ, p.zonetemp_max_occ)
            for t = 1 : o.numstages
                df_T[!,Symbol("floor$(f)_zone$(z)_T_pred_$(t)")] = [JuMP.value(zonetemp[f, z, t])]
                if  o.debugmodeflag == 1
                    new_row_input = [(minute+t), JuMP.value(zoneflow[f, z,window_index[t]]), JuMP.value(zonedischargetemp[f,z,window_index[t]]),JuMP.value(oat[t]),JuMP.value(intload[f,z,t])]
                    push!(df_input,new_row_input)
                    new_row_output = [(minute+t),JuMP.value(zonetemp[f, z, t]),JuMP.value(heatsetpoint[t]) - JuMP.value(slack[f,z,t]), JuMP.value(coolsetpoint[t]) + JuMP.value(slack[f,z,t])]
                    push!(df_output,new_row_output)
                    new_row_u = [(minute+t),JuMP.value(ahudamper[f, window_index[t]]),JuMP.value(ahusupplytemp[f, window_index[t]]), JuMP.value(zoneflow[f, z, window_index[t]])/p.zoneflow_max[f, z],
                        (JuMP.value(zonedischargetemp[f, z, window_index[t]]) - JuMP.value(ahusupplytemp[f, window_index[t]]) - 3) / (p.zonedischargetemp_max - JuMP.value(ahusupplytemp[f, window_index[t]]) - 3),
                        JuMP.value(zonedischargetemp[f, z, window_index[t]])]
                    push!(df_u, new_row_u)
                end
            end
            df_T[!,Symbol("floor$(f)_zone$(z)_T_real")] = [dfCurrentMeasurements[1, Symbol("floor$(f)_zon$(z)_TRooAir_y")]]
            @printf("floor %d, zone %d -> reheat valve opening = %.4f, constrain: [%.2f, %2f]\n", f, z, JuMP.value(zoneflow[f, z, 1])/p.zoneflow_max[f,z], p.zoneflow_min[f,z]/p.zoneflow_max[f,z], p.zoneflow_max[f,z]/p.zoneflow_max[f,z])
            if  o.debugmodeflag == 1
                CSV.write("./debuginfo/Zone$(zid)_time$(minute+1)_input.csv",df_input, append=false)
                CSV.write("./debuginfo/Zone$(zid)_time$(minute+1)_output.csv",df_output, append=false)
                CSV.write("./debuginfo/Zone$(zid)_time$(minute+1)_u.csv",df_u, append=false)
            end
        end
    end
    if o.debugmodeflag == 1
        df_C[!,Symbol("Cost1")] = [sum(JuMP.value(slack[f,z,t])^2 for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages)*o.penalty]
        df_C[!,Symbol("Cost2")] = [solverinfo["optcost"] - sum(JuMP.value(slack[f,z,t])^2 for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages)*o.penalty]
        CSV.write("./result_T_pred.csv", df_T, append=true)
        CSV.write("./result_P_pred.csv", df_P, append=true)
        CSV.write("./result_cost.csv", df_C, append=true)
    end
    global currMPCStatus = solverinfo["status"]
    if in(solverinfo["status"], status_optimalSolution)
      global currMPCflag = 0
      # store current overrides
      global currMPCStage = 1
      @printf("Using stage %d of the MPC prediction horizon.\n", currMPCStage)
      global dfCurrentSetpoints = setoverrides!(dfCurrentSetpoints, control = "MPC", stage = currMPCStage)
      @printf("<<<<<<< RE-INITIALIZE THE MODEL WITH THE LAST VALUES OF THE PREVIOUS SOLVER. >>>>>>>\n")
      JuMP.set_start_value.(JuMP.all_variables(m), JuMP.value.(JuMP.all_variables(m)))
      o.maxiter = originalMaxIter
      @printf("New max iteration number: %d.\n", o.maxiter)

    elseif in(solverinfo["status"], status_goodSolution)
      global currMPCflag = 1
      # store current overrides
      global currMPCStage = 1
      @printf("Using stage %d of the MPC prediction horizon.\n", currMPCStage)
      global dfCurrentSetpoints = setoverrides!(dfCurrentSetpoints, control = "MPC", stage = currMPCStage)
      @printf("<<<<<<< RE-INITIALIZE THE MODEL WITH THE LAST VALUES OF THE PREVIOUS SOLVER. >>>>>>>\n")
      JuMP.set_start_value.(JuMP.all_variables(m), JuMP.value.(JuMP.all_variables(m)))
      o.maxiter = originalMaxIter
      @printf("New max iteration number: %d.\n", o.maxiter)
    elseif solverinfo["status"] == status_userLim[1] # MOI.ITERATION_LIMIT
      global currMPCflag = 2
      global currMPCStage = "n/a"
      @printf("<<<<< THIS TIME: USER MAXIMUM ITERATION NUMBER REACHED. >>>>>>>\n")
      @printf("<<<<<<< RE-INITIALIZE THE MODEL WITH THE LAST VALUES OF THE PREVIOUS SOLVER. >>>>>>>\n")
      JuMP.set_start_value.(JuMP.all_variables(m), JuMP.value.(JuMP.all_variables(m)))
      # store previously computed MPC overrides
      global dfCurrentSetpoints = setoverrides!(dfCurrentSetpoints, dfPastSetpoints, minute_of_day, control = "MPC")
      o.maxiter += 100
      @printf("New max iteration number: %d.\n", o.maxiter)
    elseif in(solverinfo["status"], status_badSolution)
      global currMPCflag = 3
      global currMPCStage = "n/a"
      @printf("Optimization algorithm terminated with INFEASIBLE or UNBOUNDED solution.\n")
      @printf("<<<<<<< RE-INITIALIZE THE MODEL WITH THE LAST VALUES OF THE PREVIOUS SOLVER. >>>>>>>\n")
      JuMP.set_start_value.(JuMP.all_variables(m), JuMP.value.(JuMP.all_variables(m)))
      # store previously computed MPC overrides
      global dfCurrentSetpoints = setoverrides!(dfCurrentSetpoints, dfPastSetpoints, minute_of_day, control = "MPC")
    elseif in(solverinfo["status"], status_userLim[2:end])
      global currMPCflag = 4
      global currMPCStage = "n/a"
      @printf("<<<<<<< SOME OTHER SORT OF LIMIT HAS BEEN REACHED, that is %s. >>>>>>>>>>", string(solverinfo["status"]))
      @printf("<<<<<<< RE-INITIALIZE THE MODEL WITH THE LAST VALUES OF THE PREVIOUS SOLVER. >>>>>>>\n")
      JuMP.set_start_value.(JuMP.all_variables(m), JuMP.value.(JuMP.all_variables(m)))
      # store previously computed MPC overrides
      global dfCurrentSetpoints = setoverrides!(dfCurrentSetpoints, dfPastSetpoints, minute_of_day, control = "MPC")
    end
    mpcOptEndTime = Base.Libc.time()
    if o.debugmodeflag == 1
        df_F[!,Symbol("MPC_solution_flag")] = [currMPCflag]
        CSV.write("./result_MPC_flag.csv", df_F, append=true)
    end

    global dfPastSetpoints = deepcopy(dfCurrentSetpoints)
    # send back data (with overrides)
    u = df2dict!(u, dfCurrentSetpoints) # u
    return u
end

function initialize()
    # u: dict Defines the initial control input to be used for the next step.
    # {<input_name> : <input_value>}
    # u = Dict("oveAct_u" => 0.0,"oveAct_activate" => 1)
    # u = ctrlInitialize(case.inputs)
    global u = Dict{String, Float64}()
    for f = 1 : p.numfloors
        # damper setpoint
        u["floor$(f)_aHU_con_oveMinOAFra_activate"] = 0
        u["floor$(f)_aHU_con_oveMinOAFra_u"] = 1e-27
        # ahu supply temperatures
        u["floor$(f)_aHU_con_oveTSetSupAir_activate"] = 0
        u["floor$(f)_aHU_con_oveTSetSupAir_u"] = 35 # Celsius
        # static pressure setpoint
        u["set_ahupressure_f$(f)"] = 100
        for z = 1 : p.numzones
            u["floor$(f)_zon$(z)_oveTSetDisAir_activate"] = 0
            u["floor$(f)_zon$(z)_oveTSetDisAir_u"] = 1e-27
            u["floor$(f)_zon$(z)_oveAirFloRat_activate"] = 0
            u["floor$(f)_zon$(z)_oveAirFloRat_u"] = 1e-27
            u["floor$(f)_zon$(z)_oveHeaOut_activate"] = 0
            u["floor$(f)_zon$(z)_oveHeaOut_u"] = 1e-27
        end
    end
    # global CurrentMeasurements = Dict{String, Float64}()
    # CurrentMeasurements["TOutDryBul_y"] = 35.0
    # CurrentMeasurements["EnergyPrice"] = 1.0
    # for f = 1 : p.numfloors
    #     # ahu supply temperatures
    #     CurrentMeasurements["floor$(f)_TSupAir_y"] = 20.0
    #     CurrentMeasurements["floor$(f)_TMixAir_y"] = 25.0 # Celsius
    #     CurrentMeasurements["floor$(f)_TRetAir_y"] = 35.0 # Celsius
    #     for z = 1 : p.numzones
    #         CurrentMeasurements["floor$(f)_zon$(z)_TRooAir_y"] = 20.0 #zonetemp
    #         CurrentMeasurements["floor$(f)_zon$(z)_TSupAir_y"] = 15.0# zonedischargetemp
    #         CurrentMeasurements["floor$(f)_zon$(z)_mSupAir_y"] = 5.0 # zoneflow
    #     end
    # end
    global dfCurrentSetpoints = dict2df!(dfCurrentSetpoints, u)
    return u #, CurrentMeasurements
end

end
