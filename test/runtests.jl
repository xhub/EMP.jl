#for undef and occursin
using Compat
using Compat.Test

#examples_path = joinpath(dirname(dirname(@__FILE__)), "examples")
#for example in ["simple_mopec2.jl"]
#    include(joinpath(examples_path, example))
#end

# This doesn't appears to be in Compat?
@static if VERSION >= v"0.7.0"
    function isdef(symbol)
        return isdefined(@__MODULE__, symbol)
    end
else
    function isdef(symbol)
        return isdefined(symbol)
    end
end

# On the new macros are really robust
if VERSION >= v"1.0.0"
    include("variable.jl")
end

include("simple_mopec2.jl")
include("gnep_river_basin.jl")
include("ovf_loss_nl.jl")
include("ovf_loss_equil.jl")

# write your own tests here
