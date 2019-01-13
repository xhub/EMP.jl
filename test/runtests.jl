using EMP
using Compat.Test

#examples_path = joinpath(dirname(dirname(@__FILE__)), "examples")
#for example in ["simple_mopec2.jl"]
#    include(joinpath(examples_path, example))
#end
include("simple_mopec2.jl")
include("gnep_river_basin.jl")
include("ovf_loss_nl.jl")
include("ovf_loss_equil.jl")

# write your own tests here
