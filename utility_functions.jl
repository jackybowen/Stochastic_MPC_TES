"""
Create new directory for the current MPC run.
"""
function createdir()
    # parent directory
    parent = joinpath("./results")
    # run index
    run = 1
    # check if directory with same name exists
    while isdir(joinpath(parent, "run$run"))
        run += 1
    end
    # create directory with unique name
    dirpath = joinpath(parent, "run$run")
    mkpath(dirpath)
    return dirpath
end

"""
Receive data from client using socket connection.
"""
function getCaseInfo(url, case)
    # GET CASE INFORMATION
    # --------------------
    # Test case name
    case.name = JSON.parse(String(HTTP.get("$url/name").body))
    # Inputs available
    case.inputs = JSON.parse(String(HTTP.get("$url/inputs").body))
    # Measurements available
    case.measurements = JSON.parse(String(HTTP.get("$url/measurements").body))
    # Default simulation step
    case.step_def = JSON.parse(String(HTTP.get("$url/step").body))
    return case
end

function ctrlInitialize(inputs)
    global u = Dict{String, Float64}()
    for ind = 1 : length(inputs)
      if occursin("_u", inputs[ind])
        u[inputs[ind]] = 1e-27
      elseif occursin("_activate", inputs[ind])
        u[inputs[ind]] = 0
      end
    end
    return u
end

"""
Save current data (measurement) into a dataframe.
"""
function dict2df!(df::DataFrames.DataFrame, data::Dict)
    df = DataFrames.DataFrame([typeof(data[k]) for k in keys(data)], [Symbol(k) for k in keys(data)], 1)
    # loop over variable names in data
    for v in keys(data)
        # make sure that value is stored as float even if integer value is sent for a "Double" type sent by client
        if typeof(data[v]) == "Float64"
            # df.v = [values(data[v]) * 1.0] # store corresponding value as float
            df[1, Symbol(v)] = values(data[v]) * 1.0 # store corresponding value as float
        else
            # df.v = [values(data[v])]    # store corresponding value as  integer
            df[1, Symbol(v)] = values(data[v])    # store corresponding value as  integer
        end
    end
    return df
end

"""
Convert data frame to dict, to go from Julia setpoints to JSON u.
"""
function df2dict!(data::Dict{String, Float64}, df::DataFrames.DataFrame)
    # loop over variable names in data
    for v in keys(data)
        data[v] = df[1, Symbol(v)]
    end
    return data
end

"""
Returns the heating and cooling setpoints (comfort bounds) over the MPC
prediction horizon, starting from a given time instant.
"""
function comfortbounds(current_minute::Float64)

    h = zeros(o.numstages)   # initialize lower bounds (heating setpoints)
    c = zeros(o.numstages)   # initialize upper bounds (cooling setpoints)

    for i in 1:o.numstages
        min = current_minute + i # clock time (in minutes) for i-th stage
        hour = ceil(min/60.0)    # clock time (in hour) for i-th stage
        min_of_day = Helper.minute_of_day(min) # minute of day for i-th stage
        # println(min_of_day)
        day_of_week = Helper.day_of_week(hour) # day of the week for i-th stage

        if  o.comfortboundsflag == 1
            h[i] = true_comfortbounds[round(Int,min_of_day),:h]
            c[i] = true_comfortbounds[round(Int,min_of_day),:c]
            println(h[i]," and ",c[i])
        else
            # heating and cooling setpoints depending on time
            if day_of_week == 1.0  # Sunday
                h[i] = p.zonetemp_min_uocc
                c[i] = p.zonetemp_max_uocc
            elseif 60.0 * 5 <= min_of_day < 60.0 * 21 # occupied periods in weekdays and Saturday
                h[i] = p.zonetemp_min_occ
                c[i] = p.zonetemp_max_occ
            else #unoccupied periods in weekdays and Saturday
                h[i] = p.zonetemp_min_uocc
                c[i] = p.zonetemp_max_uocc
            end
        end
        # update hour_index with current time stamp
        hour_of_day = ceil(min_of_day / 60.0)
        global hour_index[i] =  hour_of_day

    end

    return h, c
end

"""
Predict ambient temperature over the prediction horizon.
"""
function predictambient(current_temp::Float64, minute::Float64;
                        unit::String = "Celsius")
    if o.oatpredflag == 0  # use current temp as prediction over entire horizon
        temps = current_temp * ones(o.numstages)

    elseif o.oatpredflag == 1 # use complete future information
        minute = round(Int, minute+1)   # starting row index (integer)
        minute_ind = (minute : minute + o.numstages-1)
        temps = true_outsidetemp[minute_ind.+1, :temp] # dataframe
        temps = convert(Array, temps)  # convert to Array
        temps[1] = current_temp # Use measurement instead of prediction for current step
    end

    if unit == "Kelvin"
        temps .-= 273.15  # convert current temps to Celsius from Kelvin
    end

    return temps
end

"""
Predict internal loads over the prediction horizon, given the past history.
"""
function predictloads(history::DataFrames.DataFrame,minute_of_day::Float64)

    # flag for load prediction (Int)
    flag = o.loadpredflag
    if flag == 0
        loads = zeros(p.numfloors, p.numzones, o.numstages)  # zero loads
        if p.zonetempmodelflag == 0
            for f = 1:p.numfloors, z = 1:p.numzones
                loads[f, z, :] .= p.zonetemp_q0[f,z]
            end
        end
    elseif flag == 1
        # loads computed using a moving average model
        if size(history, 1) < o.cl_MAwindow  # not enough history available
            loads = zeros(p.numfloors, p.numzones, o.numstages) # zero loads

        elseif p.zonetempmodelflag == 3
                minute = round(Int, minute_of_day)  # starting row index (integer)
                hridx = ceil(Int, minute/60.0)    # clock time (in hour) for i-th stage

                loads = zeros(p.numfloors, p.numzones, o.numstages)  # initialize
                for (f, z) in zip(1:p.numfloors, 1:p.numzones)
                    sum_error = 0.0
                    alpha = p.zonetemp_alpha3[f, z, hridx]
                    beta = p.zonetemp_beta3[f, z, hridx]
                    gamma = p.zonetemp_gamma3[f, z, hridx]
                    q3 = p.zonetemp_q3[f, z, hridx]
                    for row in 1:(size(history, 1) - 1)

                        # extract variables of interest
                        temp = history[row, Symbol("floor$(f)_zon$(z)_TRooAir_y")] # ("zonetemp_f$(f)z$z")]
                        flow = history[row, Symbol("floor$(f)_zon$(z)_mSupAir_y")] # ("zoneflow_f$(f)z$z")]
                        dischargetemp = history[row, Symbol("floor$(f)_zon$(z)_TSupAir_y")] # ("zonedischargetemp_f$(f)z$z")]
                        ambient = history[row, Symbol("TOutDryBul_y")] # ("outside_temp")]

                        # predicted temp (from model) for data in row+1 given data in row
                        pred =  alpha * temp +  beta * flow * (dischargetemp - temp) + gamma * ambient + q3
                        # actual observation for data in row+1
                        actual = history[row + 1, Symbol("floor$(f)_zon$(z)_TRooAir_y")] # ("zonetemp_f$(f)z$z")]
                        sum_error += (actual - pred)
                    end
                    loads[f, z, :] .= sum_error / o.cl_MAwindow  # average of the mispredictions
                end
        else
            loads = zeros(p.numfloors, p.numzones, o.numstages)  # initialize
            for (f, z) in zip(1:p.numfloors, 1:p.numzones)
                sum_error = 0.0
                alpha = p.zonetemp_alpha[f, z]
                beta = p.zonetemp_beta[f, z]
                gamma = p.zonetemp_gamma[f, z]
                for row in 1:(size(history, 1) - 1)

                    # extract variables of interest
                    temp = history[row, Symbol("floor$(f)_zon$(z)_TRooAir_y")] # ("zonetemp_f$(f)z$z")]
                    flow = history[row, Symbol("floor$(f)_zon$(z)_mSupAir_y")] # ("zoneflow_f$(f)z$z")]
                    dischargetemp = history[row, Symbol("floor$(f)_zon$(z)_TSupAir_y")] # ("zonedischargetemp_f$(f)z$z")]
                    ambient = history[row, Symbol("TOutDryBul_y")] # ("outside_temp")]

                    # predicted temp (from model) for data in row+1 given data in row
                    pred =  alpha * temp +  beta * flow * (dischargetemp - temp) + gamma * ambient
                    # actual observation for data in row+1
                    actual = history[row + 1, Symbol("floor$(f)_zon$(z)_TRooAir_y")] # ("zonetemp_f$(f)z$z")]
                    sum_error += (actual - pred)
                end
                loads[f, z, :] .= sum_error / o.cl_MAwindow  # average of the mispredictions
            end
        end
    elseif flag == 2
        loads = zeros(p.numfloors, p.numzones, o.numstages+1)  # initialize
        # println("Outside air temperature: ",history[end, Symbol("TOutDryBul_y")])
        # println("Hour:",hour," Minute:",minute_of_day)
        for f=1:p.numfloors, z=1:p.numzones, t= 0:(o.numstages-1)
            minute = round(Int, minute_of_day) + t   # starting row index (integer)
            hour = ceil(Int, minute/60.0)    # clock time (in hour) for i-th stage
            if t == 0
               loads[f,z,1] = p.zonetemp_b[f,z,hour+2] + p.zonetemp_b[f,z,2]*history[end, Symbol("TOutDryBul_y")] + p.zonetemp_b[f,z,1]
            else
            # println(hour)
               loads[f,z,t+1] = p.zonetemp_b[f,z,hour+2] + p.zonetemp_b[f,z,2]*true_outsidetemp[minute,:temp] + p.zonetemp_b[f,z,1]
            end
        end
    end

    return loads
end

"""
Set appropriate setpoints from the MPC model.
"""
function setoverrides!(df::DataFrames.DataFrame;
                        control = "MPC",
                        stage::Int64 = 1,
                        default::Float64 = 1e-27,
                        unit::String = "Celsius")
    if control == "MPC"
        for f = 1:p.numfloors
                # damper setpoint
                df[1, Symbol("floor$(f)_aHU_con_oveMinOAFra_activate")] = 1
                df[1, Symbol("floor$(f)_aHU_con_oveMinOAFra_u")] = JuMP.value(ahudamper[f, stage])

                # ahu supply temperatures
                df[1, Symbol("floor$(f)_aHU_con_oveTSetSupAir_activate")] = 1
                if unit == "Kelvin"
                    df[1, Symbol("floor$(f)_aHU_con_oveTSetSupAir_u")] = JuMP.value(ahusupplytemp[f, stage]) + 273.15
                elseif unit == "Celsius"
                    df[1, Symbol("floor$(f)_aHU_con_oveTSetSupAir_u")] = JuMP.value(ahusupplytemp[f, stage])
                end

                # static pressure setpoint
                mflow = JuMP.value(sum_zoneflows[f, 1])
                # if in(Symbol("set_ahupressure_f$(f)"), names(df))
                #     df[1, Symbol("set_ahupressure_f$f")] = staticpressure(mflow)
                # else
                #     insertcols!(df, size(df, 2) + 1, Symbol("set_ahupressure_f$f") => staticpressure(mflow))
                # end
                df[1, Symbol("set_ahupressure_f$f")] = staticpressure(mflow, f)

                ## zone-level setpoints
                for z = 1:p.numzones
                    # zone flows
                    df[1, Symbol("floor$(f)_zon$(z)_oveAirFloRat_activate")] = 1
                    df[1, Symbol("floor$(f)_zon$(z)_oveAirFloRat_u")] = JuMP.value(zoneflow[f, z, stage])/p.zoneflow_max[f, z]

                    df[1, Symbol("floor$(f)_zon$(z)_oveHeaOut_activate")] = 1
                    df[1, Symbol("floor$(f)_zon$(z)_oveHeaOut_u")] = (JuMP.value(zonedischargetemp[f, z, stage]) - JuMP.value(ahusupplytemp[f, stage])-3) / (p.zonedischargetemp_max - JuMP.value(ahusupplytemp[f, stage])-3)

                    # discharge temperatures - this is just for data saving purpose
                    # and it comes in Celsius from the MPC algorithm
                    df[1, Symbol("floor$(f)_zon$(z)_oveTSetDisAir_activate")] = 1
                    df[1, Symbol("floor$(f)_zon$(z)_oveTSetDisAir_u")] = JuMP.value(zonedischargetemp[f, z, stage])
                end
        end
    elseif control == "DEFAULT"
        for f = 1:p.numfloors
            # damper setpoint
            df[1, Symbol("floor$(f)_aHU_con_oveMinOAFra_activate")] = 0
            df[1, Symbol("floor$(f)_aHU_con_oveMinOAFra_u")] = default
            # ahu supply temperatures
            df[1, Symbol("floor$(f)_aHU_con_oveTSetSupAir_activate")] = 0
            df[1, Symbol("floor$(f)_aHU_con_oveTSetSupAir_u")] = default

            # static pressure setpoint
            if in(Symbol("set_ahupressure_f$(f)"), names(df))
                df[1, Symbol("set_ahupressure_f$f")] = default
            else
                insertcols!(df, size(df, 2) + 1, Symbol("set_ahupressure_f$f") => default)
            end
            # df[1, Symbol("set_ahupressure_f$f")] = default

            ## zone-level setpoints
            for z = 1:p.numzones
                # zone flows
                df[1, Symbol("floor$(f)_zon$(z)_oveAirFloRat_activate")] = 0
                df[1, Symbol("floor$(f)_zon$(z)_oveAirFloRat_u")] = default

                df[1, Symbol("floor$(f)_zon$(z)_oveHeaOut_activate")] = 0
                df[1, Symbol("floor$(f)_zon$(z)_oveHeaOut_u")] = default

                # discharge temperatures - this is just for data saving purpose
                df[1, Symbol("floor$(f)_zon$(z)_oveTSetDisAir_activate")] = 0
                df[1, Symbol("floor$(f)_zon$(z)_oveTSetDisAir_u")] = default
            end
        end
    else
        nothing  # error message
    end
    return df
end

"""
Copy setpoint values in a target dataframe from a source dataframe.
"""
function setoverrides!(target::DataFrames.DataFrame,
                       source::DataFrames.DataFrame,
                       minute_of_day::Float64;
                       control = "MPC")
    if control == "MPC"
        target = source
    else
        nothing # display error message
    end

    return target
end

"""
Update MPC parameters using predicted outside temp, predicted internal loads
and comfort bounds.
"""
function updatemodelparams(df::DataFrames.DataFrame, params::Dict)
    # extract relevant parameters
    h = params["heat_sp"]
    c = params["cool_sp"]
    pred_oat = params["pred_oat"]
    pred_loads = params["pred_loads"]
################################################################################
    JuMP.set_value(price, df[1, Symbol("EnergyPrice")])
################################################################################

    for tInd in 1:o.numstages
      # update outside-air temperature parameters
      JuMP.set_value(oat[tInd], pred_oat[tInd])
      # update heating and cooling setpoints
      JuMP.set_value(heatsetpoint[tInd], h[tInd])
      JuMP.set_value(coolsetpoint[tInd], c[tInd])
      # Tsa - Tma is unconstrained for all stages
      JuMP.set_value(bound[tInd], 1e2)
    end

    # update internal-load parameters
    for fInd = 1:p.numfloors, zInd = 1:p.numzones, tInd = 1:o.numstages
        JuMP.set_value(intload[fInd, zInd, tInd], pred_loads[fInd, zInd, tInd])
    end
    if p.zonetempmodelflag == 0
        # update initial zone temperature
        for fInd = 1:p.numfloors, zInd = 1:p.numzones
            JuMP.set_value(zonetemp_init[fInd, zInd], df[1, Symbol("floor$(fInd)_zon$(zInd)_TRooAir_y")]) #  "zonetemp_f$(fInd)z$(zInd)")])
        end
    end
    if p.adpflag == 1
        datafile = "./a_parameters_adaptive.csv"
        for f in 1:p.numfloors, z in 1:p.numzones
            p.zonetemp_alpha[f, z] = CSV.read(datafile,DataFrame;header=0)[!,Symbol(string("Column",1))][(f-1)*(p.numzones) + z]
            p.zonetemp_beta[f, z] = CSV.read(datafile,DataFrame;header=0)[!,Symbol(string("Column",2))][(f-1)*(p.numzones) + z]
            p.zonetemp_gamma[f, z] = CSV.read(datafile,DataFrame;header=0)[!,Symbol(string("Column",3))][(f-1)*(p.numzones) + z]
            p.zonetemp_q0[f, z] = CSV.read(datafile,DataFrame;header=0)[!,Symbol(string("Column",4))][(f-1)*(p.numzones) + z]
        end
    end
    # first set zone flows and ahudampers to correct bounds (if it was currently in state 1)
    for fInd = 1:p.numfloors, zInd = 1:p.numzones, tInd = 1:o.numstages
        hInd = ceil(Int64, tInd / o.controlwindow);
        if  5 <= hour_index[tInd] < 21
            JuMP.set_lower_bound(zoneflow[fInd, zInd, hInd], p.zoneflow_min[fInd, zInd])
            JuMP.set_lower_bound(ahudamper[fInd, hInd],p.damper_min[fInd])
        else
            JuMP.set_lower_bound(zoneflow[fInd, zInd, hInd], 0.0)
            JuMP.set_lower_bound(ahudamper[fInd, hInd], 0.0)
        end
        JuMP.set_upper_bound(zoneflow[fInd, zInd, hInd], p.zoneflow_max[fInd, zInd])
    end
end

function updatemodelinit(history::DataFrames.DataFrame,minute_of_day::Float64)
    # update initial zone temperature
    minute = round(Int, minute_of_day) - 1   # starting row index (integer)
    hour = ceil(Int, minute/60.0)    # clock time (in hour) for i-th stage
    if size(history,1) <= 1
       Tout = history[end, Symbol("TOutDryBul_y")]
    else
       Tout = history[end-1, Symbol("TOutDryBul_y")]
    end

    for fInd = 1:p.numfloors, zInd = 1:p.numzones
        exist_load = p.zonetemp_b[fInd,zInd,hour+2] + p.zonetemp_b[fInd,zInd,2]*Tout + p.zonetemp_b[fInd,zInd,1]
        JuMP.set_value(zonetemp_init[fInd, zInd], history[end, Symbol("floor$(fInd)_zon$(zInd)_TRooAir_y")]-exist_load) #  "zonetemp_f$(fInd)z$(zInd)")])
        # println("floor",fInd," zone",zInd," zone flow: ",df[1, Symbol("floor$(fInd)_zon$(zInd)_mSupAir_y")])
        JuMP.set_value(zonedischargetemp_init[fInd, zInd], history[end, Symbol("floor$(fInd)_zon$(zInd)_TSupAir_y")]) #  "zonetemp_f$(fInd)z$(zInd)")])
        JuMP.set_value(zoneflow_init[fInd, zInd], max(0,history[end, Symbol("floor$(fInd)_zon$(zInd)_mSupAir_y")])) #  "zonetemp_f$(fInd)z$(zInd)")])
    end
end

"""
Update table of input data over a rolling horizon.
"""
function updatehistory!(history::DataFrames.DataFrame, current::DataFrames.DataFrame)
    history = vcat(history, current)     # append current measurements at the end of history
    if size(history, 1) > o.cl_MAwindow  # check if number of rows > length of rolling window
        history = history[2:end, :]      # remove first row
    end
    return history
end

"""
Solve MPC model.
"""
function solvemodel()
    # message
    println("Solving MPC model ... ")

    # solve mpc
    _, solTime, solAllocatedBytes, garbageCollectorTime, solMemAllocs = @timed JuMP.optimize!(m)
    status = JuMP.termination_status(m)
    objValue = JuMP.objective_value(m)

    # solution info
    solverinfo = Dict("optcost"   => objValue,
                      "status"    => status,
                      "soltime"   => solTime)
    return solverinfo
end

"""
Determine static pressure setpoints using convex hull information.
"""
function staticpressure(mflow::Float64, fInd::Int64)
    # # massflow breakpoints
    # m0, m1, m2, m3, m4, m5 = sum(p.zoneflow_min), 5.81, 16.68, 17.05, 17.61, sum(p.zoneflow_max)
    # # pressure breakpoints
    # p0, p1, p2, p3, p4, p5 = 24.88, 24.88, 121.51, 128.68, 160.88, 160.88
    f = fInd # Index of floor level
    # massflow breakpoints
    m0, m1, m2, m3, m4, m5 = p.massflow_sample_0[f], p.massflow_sample_1[f], p.massflow_sample_2[f], p.massflow_sample_3[f], p.massflow_sample_4[f], p.massflow_sample_5[f]
    # pressure breakpoints
    p0, p1, p2, p3, p4, p5 = p.pressure_sample_0[f], p.pressure_sample_1[f], p.pressure_sample_2[f], p.pressure_sample_3[f], p.pressure_sample_4[f], p.pressure_sample_5[f]

    # piecewise-linear pressure function
    if  m0 <= mflow <= m1
        pressure = p0 + (p1 - p0) * (mflow - m0) / (m1 - m0)
    elseif m1 < mflow <= m2
        pressure = p1 + (p2 - p1) * (mflow - m1) / (m2 - m1)
    elseif m2 < mflow <= m3
        pressure = p2 + (p3 - p2) * (mflow - m2) / (m3 - m2)
    elseif m3 < mflow <= m4
        pressure = p3 + (p4 - p3) * (mflow - m3) / (m4 - m3)
    else
        pressure = p4 + (p5 - p4) * (mflow - m4) / (m5 - m4)
    end
    return pressure
end
