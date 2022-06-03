module LoadParams
# include standard packages
using ConfParser
using CSV

"""
Function to parse strings from config files to integer values.
"""
function parse_int(conf::ConfParser.ConfParse, section::AbstractString, variable::AbstractString)
    x = ConfParser.retrieve(conf, section, variable)
    return round(Int64, parse(Float64, x))
end

"""
Function to parse strings from config files to floating-point values.
"""
function parse_float(conf::ConfParser.ConfParse, section::AbstractString, variable::AbstractString)
    x = ConfParser.retrieve(conf, section, variable)
    return parse(Float64, x)
end

"""
Function to parse strings from config files to boolean values.
"""
function parse_bool(conf::ConfParser.ConfParse, section::AbstractString, variable::AbstractString)
    x = ConfParser.retrieve(conf, section, variable)
    return parse(Bool, x)
end

# define user-defined data types to store model parameters
"""
Data type to store building-specific parameters.
"""
mutable struct BldgParams
    numzones::Int
    numfloors::Int
    eff_cool::Float64
    eff_heat::Float64
    specheat::Float64
    modelflag::Int
    season_month::Int

    zonetemp_min_occ::Float64              # heating setpoint in occupied period
    zonetemp_max_occ::Float64              # cooling setpoint in occupied period
    zonetemp_min_uocc::Float64             # heating setpoint in unoccupied period
    zonetemp_max_uocc::Float64             # cooling setpoint in unoccupied period

    zoneflow_min::Array{Float64,2}         # of size (numfloors, numzones)
    zoneflow_max::Array{Float64,2}         # of size (numfloors, numzones)

    zonedischargetemp_max::Float64         # max zone discharge temp max
    ahusupplytemp_min::Float64
    ahusupplytemp_max::Float64
    ahusupplytemp_max_dev::Float64
    ahuflow_min::Vector{Float64}           # of length numfloors
    ahuflow_max::Vector{Float64}           # of length numfloors

    fan_params::Array{Float64,2}            # of length 4
    # pressure_min::Vector{Float64}
    # pressure_max::Vector{Float64}
    damper_min::Vector{Float64}
    damper_max::Vector{Float64}

    # inner constructor to initilize an object of type Params
    function BldgParams()
        # read config file for optimization and model parameters
        conf = ConfParse("Configs/bldg_params.ini")
        parse_conf!(conf)

        # define new instance
        obj = new()

        ## populate the fields of the instance
        # building-specific params
        section = "building_params"
        obj.numzones = parse_int(conf, section, "numzones")
        obj.numfloors = parse_int(conf, section, "numfloors")
        obj.eff_cool = parse_float(conf, section, "eff_cool")
        obj.eff_heat = parse_float(conf, section, "eff_heat")
        obj.specheat = parse_float(conf, section, "specheat")
        obj.modelflag = parse_int(conf, section, "modelflag")
        obj.season_month = parse_int(conf, section, "season_month")

        # comfort-level params
        section = "comfort_params"
        obj.zonetemp_min_occ = parse_float(conf, section, "zonetemp_min_occ")
        obj.zonetemp_max_occ = parse_float(conf, section, "zonetemp_max_occ")
        obj.zonetemp_min_uocc = parse_float(conf, section, "zonetemp_min_uocc")
        obj.zonetemp_max_uocc = parse_float(conf, section, "zonetemp_max_uocc")

        # zone flow parameters
        section = "zoneflow_params"
        obj.zoneflow_min = ones(obj.numfloors, obj.numzones)
        obj.zoneflow_max = ones(obj.numfloors, obj.numzones)
        for f in 1:obj.numfloors, z in 1:obj.numzones
            obj.zoneflow_min[f, z] = parse_float(conf, section, "floor$(f)_zone$(z)_flow_min")
            obj.zoneflow_max[f, z] = parse_float(conf, section, "floor$(f)_zone$(z)_flow_max")
        end

        # zone discharge temp parameter
        section = "dischargetemp"
        obj.zonedischargetemp_max = parse_float(conf, section, "zonedischargetemp_max")

        # AHU level parameters
        section = "ahu_params"
        obj.ahusupplytemp_min = parse_float(conf, section, "ahusupplytemp_min")
        obj.ahusupplytemp_max = parse_float(conf, section, "ahusupplytemp_max")
        obj.ahusupplytemp_max_dev = parse_float(conf, section, "ahusupplytemp_max_dev")
        obj.ahuflow_min = ones(obj.numfloors)
        obj.ahuflow_max = ones(obj.numfloors)
        for f in 1:obj.numfloors
            obj.ahuflow_min[f] = parse_float(conf, section, "floor$(f)_ahuflow_min")
            obj.ahuflow_max[f] = parse_float(conf, section, "floor$(f)_ahuflow_max")
        end

        # fan parameters
        section = "fan_params"
        obj.fan_params = ones(obj.numfloors, 4)
        for f in 1:obj.numfloors
            for i in 1:4
                obj.fan_params[f, i] = parse_float(conf, section, "floor$(f)_fan_param$(i)")
            end
        end

        return obj
    end
end

"""
Data type to store optimization model and solver parameters.
"""
mutable struct MPCParams
    numstages::Int
    controlwindow::Int
    numwindows::Int
    weight_fan::Float64
    weight_heat::Float64
    weight_cool::Float64
    weight_norm::Float64
    penalty::Float64
    objectivemodelflag::Int
    oatpredflag::Int
    comfortboundsflag::Int
    loadpredflag::Int
    debugmodeflag::Int
    solverflag::Int
    maxiter::Int
    # cl_baseline::Int
    # cl_startday::Float64
    # cl_numdays::Float64
    # cl_nompcdays::Float64
    # cl_nosolvewindow::Int
    cl_MAwindow::Int
    # cl_rate_supplytemp::Float64
    # cl_minPerSample::Float64
    # mpcMovingBlockImpl::Bool

    # inner constructor
    function MPCParams()
        # read config file for optimization and model parameters
        conf = ConfParse("Configs/opt_params.ini")
        parse_conf!(conf)

        # define new instance
        obj = new()

        ## populate the fields of the instance
        # model params
        section = "model_params"
        obj.numstages = parse_int(conf, section, "numstages")
        obj.controlwindow = parse_int(conf, section, "controlwindow")
        obj.numwindows = obj.numstages % obj.controlwindow == 0 ?
                        round(Int64, obj.numstages / obj.controlwindow) :
                        error("Number of stages not divisible by control window; choose new values")
        obj.weight_fan = parse_float(conf, section, "weight_fan")
        obj.weight_heat = parse_float(conf, section, "weight_heat")
        obj.weight_cool = parse_float(conf, section, "weight_cool")
        obj.weight_norm = parse_float(conf, section, "weight_norm")
        obj.penalty = parse_float(conf, section, "penalty")

        # flags
        section = "flags"
        obj.objectivemodelflag = parse_int(conf, section, "objectivemodelflag")
        obj.oatpredflag = parse_int(conf, section, "oatpredflag")
        obj.loadpredflag = parse_int(conf, section, "loadpredflag")
        obj.comfortboundsflag = parse_int(conf, section, "comfortboundsflag")
        obj.debugmodeflag = parse_int(conf, section, "debugmodeflag")

        # solver params
        section = "solver_params"
        obj.solverflag = parse_int(conf, section, "solverflag")
        obj.maxiter = parse_int(conf, section, "maxiter")

        # # closedloop params
        obj.cl_MAwindow = parse_int(conf, "other_params", "MAwindow")
        return obj
    end
end
end
