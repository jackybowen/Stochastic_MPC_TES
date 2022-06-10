# include standard packages
using DataFrames
using CSV
using JuMP
using MathOptInterface
using Ipopt
## Load parameters from configs folder
include("LoadParams.jl")
# initialize instances storing the required parameters
o = LoadParams.MPCParams()        # instance of model-and-solver-specific parameters
p = LoadParams.BldgParams()       # instance of building-specific parameters

## Load utility functions for Optimization solver
include("utility_functions.jl")
## Build the initial MPC model in JuMP. Some parameters are initialized
# using default values which are dynamically updated during closed-loop runs

# define Optimization model
m = JuMP.Model(with_optimizer(Ipopt.Optimizer, max_iter = o.maxiter, print_level = 3))
set_time_limit_sec(m, 300.0)
# map decision stages to control windows in MPC
window_index = Dict(t => round(Int, ceil(t/o.controlwindow)) for t = 1:o.numstages)
hour_index = Dict(t => round(Int, 1) for t = 1:o.numstages)
trueoat = DataFrames.DataFrame()
## Add Model constraints and objective terms per components
if o.VBModel_flag == 1
    include("ModelVB.jl")
    ModelVB!(m, o, p)
end
if o.TESModel_flag == 1
    include("ModelTES.jl")
    ModelTES!(m, o, p)
end
if o.BESModel_flag == 1
    include("ModelBES.jl")
    ModelBES!(m, o, p)
end
# if  o.oatpredflag == 1
#     if  p.season_month == 4
#         datafile = "daily_oat_april.csv"
#     elseif p.season_month == 8
#         datafile = "daily_oat_august.csv"
#     end
#     trueoat = CSV.read(datafile, DataFrame);
# end
# global true_comfortbounds = DataFrames.DataFrame()
# if  o.comfortboundsflag == 1
#     true_comfortbounds = CSV.read("Configs/daily_comfortbounds.csv", DataFrame);
# end
## Objective function
@NLobjective(m, Min, 60.0 * price * sum(o.weight_heat * heatingcapacity[f, t] + o.weight_cool *chillercapacity[f,t] +
                o.penalty * sum(slack[f,z,t]^2 for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages) for f = 1:p.numfloors, z = 1:p.numzones, t = 1:o.numstages))
