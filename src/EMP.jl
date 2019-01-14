__precompile__()

module EMP

using JuMP
using ReSHOP

#for 0.6

if VERSION >= v"0.7"
    # quickfix for Nullable
    using Nullables
else
    using Compat.Iterators.flatten
    using Compat
end

export @variableMP, @objectiveMP, @NLobjectiveMP, @constraintMP, @NLconstraintMP, @vipair, @NLvipair, solveEMP, _solveEMP, MathPrgm, EquilibriumProblem, BilevelProblem, addvar!, addequ!, addovf!, getsolution, status, get_solve_result, get_solve_result_num, get_model_result, get_model_result_num, get_solve_message, get_solve_exitcode, solve, getobjval

" Mathematical Programm representation "
mutable struct MathPrgm
    emp::Any
    vars::Vector{Int}
    equs::Vector{Int}
    matching::Dict{Int,Int}
    objequ::Int
    objvar::Int
    sense::Symbol
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    solverobj::Ptr{}
end

" Optimal Value Function (OVF) representation "
mutable struct OVF
    vidx::Int
    args::Vector{Int}
    name::String
    params::Dict{String,Any}
    # No idea how to enforce params to be a Dict{String,Any}
    function OVF(v::JuMP.Variable, args::Vector{JuMP.Variable}, name::String, params::Dict)
        args_idx = [arg.col for arg in args]
        new(v.col, args_idx, name, params)
    end
end


" EMP master object "
mutable struct Model
    model_ds
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    root::Nullable{Union{MathPrgm, Vector{MathPrgm},Vector{Vector{MathPrgm}}}}
    equs::Vector{Tuple{Symbol, Int}}
    ovfs::Vector{OVF}
end

"""
    Model([modeling_pkg = "JuMP", solver = "jams"])

Create an EMP master object, using the modeling package given as argument to store the variables and equations.
The solver argument is used to construct the ReSHOPSolver object.
"""
function Model(; modeling_pkg = "JuMP", solver = "jams")
    if modeling_pkg == "JuMP"
        solver_reshop = ReSHOP.ReSHOPSolver(solver)
        model_ds = JuMP.Model(solver=solver_reshop)
        emp = Model(model_ds, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), Nullable{Union{MathPrgm,Vector{MathPrgm},Vector{Vector{MathPrgm}}}}(), Vector{Tuple{Symbol, Int}}(), Vector{OVF}())
    elseif modeling_pkg == "Convex"
        error("Convex.jl is WIP")
    else
        error("No valid modeling_ds value passed. It should be either ``JuMP'' or ``Convex''")
    end

    setemp!(model_ds, emp)

    return emp
end

"""
    Model(model_ds)

Create an EMP master object and use the argument as the model data storage object
"""
function Model(model_ds)
    emp = Model(model_ds, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), Nullable{Union{MathPrgm,Vector{MathPrgm},Vector{Vector{MathPrgm}}}}(), Vector{Tuple{Symbol, Int}}(), Vector{OVF}())

    setemp!(model_ds, emp)

    return emp
end

"""
    MathPrgm(m::EMP.Model)

Create a Mathematical Programm in the EMP master object
"""
function MathPrgm(m::EMP.Model)
    mp = MathPrgm(m, Vector{Int}(), Vector{Int}(), Dict{Int,Int}(), -1, -1, :undef, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(),C_NULL)
    push!(m.mps, mp)
    return mp
end

" Create an Equilibrium object (not implemented) "
function Equilibrium(emp::EMP.Model, mps::Vector{MathPrgm})
    push!(emp.equils, mps)
end

"""
    EquilibriumProblem(emp::EMP.Model, mps::Vector{MathPrgm})

Define the EMP problem as an equilibrium problem
"""
function EquilibriumProblem(emp::EMP.Model, mps::Vector{MathPrgm})
    emp.root = mps
end

function BilevelProblem(emp::EMP.Model, upper::MathPrgm, lower::MathPrgm)
    emp.root = upper
    push!(upper.mps, lower)
end

include("helpers.jl")

"""
    addvar!(mp, var)

Add a variable to a Mathematical Programm
"""
function addvar!(mp::MathPrgm, var::JuMP.Variable)
    push!(mp.vars, var.col)
end

function addvar!(mp::MathPrgm, var::Array{JuMP.Variable})
    map(x -> addvar!(mp, x), flatten(var))
end

function addvar!(mp::MathPrgm, var::JuMP.JuMPArray)
    map(x-> addvar!(mp, x), values(var))
end

"""
    addequ!(mp, eqn)

Add an equation to a Mathematical Programm. The argument eqn is a single or list of JuMP.ConstraintRef, not an expression
"""
function addequ!(mp::MathPrgm, eqn::JuMP.ConstraintRef)
    gidx = reg_equ(mp.emp, eqn)
    push!(mp.equs, gidx)
    return gidx
end

function addequ!(mp::MathPrgm, eqn::Vector{ConstraintRef{M,C}}) where {M <: JuMP.AbstractModel, C <: JuMP.AbstractConstraint}
    map(x -> addequ!(mp, x), eqn)
end

function addequ!(mp::MathPrgm, eqn::Vector{JuMP.ConstraintRef})
    map(x -> addequ!(mp, x), eqn)
end

function addequ!(mp::MathPrgm, eqn)
    if isa(eqn, Array)
        map(x -> addequ!(mp, x), eqn)
    else
        throw("Unknown equation type $(typeof(eqn)) for argument $eqn")
    end

end

function addequ!(mp::MathPrgm, eqn::JuMP.JuMPArray{JuMP.ConstraintRef, 1, Tuple{Int64}})
    map(x-> addequ!(mp, x), values(eqn))
end

@static if VERSION < v"0.7"
    include("macros_v0.6.jl")
else
    include("macros.jl")
end




"""
     addovf!(emp, v, args, name, params)

Tag a variable as a supported OVF

# Arguments
- `emp::EMP.Model`: the EMP object
- `v::JuMP.Variable`: the variable
- `args::Vector{JuMP.Variable}`: the arguments for the OVF problem
- `name::String`: the name of the OVF
- `params::Dict`: the parameters for defining this OVF
"""
function addovf!(emp::EMP.Model, v::JuMP.Variable, args::Vector{JuMP.Variable}, name::String, params::Dict)
    push!(emp.ovfs, OVF(v, args, name, params))
end

"""
    solveEMP(emp)

Solve an EMP model
"""
function solveEMP(emp::EMP.Model)
    JuMP.solve(emp.model_ds)
end

#########################################################################################################
# End of public API
#########################################################################################################

include("emptree.jl")

function _solveEMP(emp::EMP.Model)

    # Steps:
    # 1. Define the full models and agents
    # 2. Solve
    # 3. Report values

    reshop_model = getReSHOPModel(emp)

    ReSHOP.make_var_index!(reshop_model)
#    ReSHOP.make_con_index!(reshop_model)

    ########
    #
    # We have to mirror the behavior of ReSHOP, where the nonlinear equation are put before
    # the linear ones.
    #
    # TODO(xhub) fix this mess ...
    #
    #####

    if (length(emp.mps) > 0)
        make_equ_index!(emp.equs, emp.model_ds.internalModel)
    else
        ReSHOP.make_con_index!(reshop_model)
    end

    ctx = ReSHOP.create_reshop_ctx(reshop_model)
    reshop_model.reshop_ctx = ctx

    reshop_emp = ReSHOP.emp_create(ctx)

    # Define the EMP tree
    emptree(emp, reshop_emp, ctx)

    for ovf in emp.ovfs
        ReSHOP.reshop_ovf(reshop_emp, ovf)
    end

    if (length(emp.mps) > 0)
        ReSHOP.reshop_set_modeltype(ctx, ReSHOP.emp)
    else
        ReSHOP.reshop_set_modeltype(reshop_model)
    end

    setvarnames(ctx, emp.model_ds)

    gams_ctx, gams_dir = ReSHOP.reshop_setup_gams()
    reshop_model.reshop_ctx_dest = gams_ctx
    reshop_model.gams_dir = gams_dir

    reshop_model.solve_exitcode = ReSHOP.reshop_solve(ctx, gams_ctx, reshop_model.solver_name, reshop_emp)

    ###############################################################################################
    #
    #  Report values and codes
    #
    ###############################################################################################

    if reshop_model.solve_exitcode == 0
        ReSHOP.emp_report_values(reshop_emp, ctx)
        ReSHOP.report_results_common(reshop_model)
        _setvalues(emp.model_ds, getsolution(emp))
    else
        println("DEBUG: solver failed with status $(reshop_model.solve_exitcode)")
        reshop_model.status = :Error
        reshop_model.solution = fill(NaN,reshop_model.nvar)
        reshop_model.solve_result = "failure"
        reshop_model.solve_result_num = 999
    end

    #    TODO(xhub) move that elsewhere
#    ReSHOP.emp_delete(reshop_emp)
end

###################################################################################################
#
#   Accessor and such
#
###################################################################################################

getsolution(emp::EMP.Model) = copy(getReSHOPModel(emp).solution)
getsolvetime(emp::EMP.Model) = getReSHOPModel(emp).solve_time
status(emp::EMP.Model) = status(getReSHOPModel(emp))

function getobjval(mp::MathPrgm)
    reshop_model = getReSHOPModel(mp.emp)
    objvar = ReSHOP.emp_mp_getobjvar(mp.solverobj)
    ReSHOP.ctx_getvarval(reshop_model.reshop_ctx, objvar)
end

# Access to solve results
get_solve_result(emp::EMP.Model) = getReSHOPModel(emp).solve_result
get_solve_result_num(emp::EMP.Model) = getReSHOPModel(emp).solve_result_num
get_model_result(emp::EMP.Model) = getReSHOPModel(emp).model_result
get_model_result_num(emp::EMP.Model) = getReSHOPModel(emp).model_result_num
get_solve_message(emp::EMP.Model) = getReSHOPModel(emp).solve_message
get_solve_exitcode(emp::EMP.Model) = getReSHOPModel(emp).solve_exitcode

function Base.show(io::IO, m::EMP.Model)
    print(io, "EMP data structure")
end

end
