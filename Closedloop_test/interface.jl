cd(dirname(@__FILE__))
include("opt_main_det_test.jl")
# include("opt_main_sto_test.jl")

using JSON
using CSV
using .Controller
################### generate dict; no be need in our tool ######################

GUI_input = JSON.parsefile(ARGS[1])  # parse and transform data
println(GUI_input)
########################### call optimization ##################################

json_file_result = ARGS[2]

result = Controller.run(GUI_input)

open(json_file_result, "w") do f
    JSON.print(f, result, 4)
end
