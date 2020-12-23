__precompile__()

module EMP

using JuMP
using ReSHOP
using MathOptInterface
MOI = MathOptInterface

using Compat
using Base.Iterators

export @vipair, @NLvipair, EMPmaster, solveEMP, _solveEMP, MathPrgm, EquilibriumProblem, BilevelProblem, addovf!, getsolution, status, termination_status, primal_status, dual_status, get_model_result_num, solve, getobjval


include("JuMPExtension.jl")

" Optimal Value Function (OVF) representation "
mutable struct OVF
    vidx::Int
    args::Vector{Int}
    name::String
    params::Dict{String,Any}
    # No idea how to enforce params to be a Dict{String,Any}
    function OVF(v::JuMP.VariableRef, args::Vector{JuMP.VariableRef}, name::String, params::Dict)
        args_idx = [arg.index.value for arg in args]
        new(v.index.value, args_idx, name, params)
    end
end


" EMP master object "
mutable struct EMPmaster
    backend
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    root::Union{Nothing,MathPrgm,Vector{MathPrgm},Vector{Vector{MathPrgm}}}
    ovfs::Vector{OVF}
    NLoffset::Int
end

"""
    Model([modeling_pkg = "JuMP", solver = "jams"])

Create an EMP master object, using the modeling package given as argument to store the variables and equations.
The solver argument is used to construct the ReSHOPSolver object.
"""
function EMPmaster(; modeling_pkg = "JuMP", solver = "jams")
    if modeling_pkg == "JuMP"
        backend = direct_model(ReSHOP.Optimizer(;solver=solver))
        emp = EMPmaster(backend)
    elseif modeling_pkg == "Convex"
        error("Convex.jl is WIP")
    else
        error("No valid modeling_ds value passed. It should be either ``JuMP'' or ``Convex''")
    end

    return emp
end

"""
    Model(backend)

Create an EMP master object and use the argument as the model data storage object
"""
function EMPmaster(backend)
    if !(backend.moi_backend isa ReSHOP.Optimizer)
        error("The JuMP model must be direct solver, use JuMP.direct_model")
    end
    emp = EMPmaster(backend, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), nothing, Vector{OVF}(), -1)

    return emp
end

"""
    MathPrgm(m::EMPmaster)

Create a Mathematical Programm in the EMP master object
"""

" Create an Equilibrium object (not implemented) "
function Equilibrium(emp::EMPmaster, mps::Vector{MathPrgm})
    push!(emp.equils, mps)
end

"""
    EquilibriumProblem(emp::EMPmaster, mps::Vector{MathPrgm})

Define the EMP problem as an equilibrium problem
"""
function EquilibriumProblem(emp::EMPmaster, mps::Vector{MathPrgm})
    emp.root = mps
end

function BilevelProblem(emp::EMPmaster, upper::MathPrgm, lower::MathPrgm)
    emp.root = upper
    push!(upper.mps, lower)
end

include("helpers.jl")

#"""
#    addvar!(mp, var)
#
#Add a variable to a Mathematical Programm
#"""
#function addvar!(mp::MathPrgm, var::JuMP.VariableRef)
#    push!(mp.vars, var.index.value)
#end
#
#function addvar!(mp::MathPrgm, var::Array{JuMP.VariableRef})
#    map(x -> addvar!(mp, x), flatten(var))
#end
#
#function addvar!(mp::MathPrgm, var::JuMP.JuMPArray)
#    map(x-> addvar!(mp, x), values(var))
#end
#
#function addvar!(mp::MathPrgm, var::JuMP.JuMPDict)
#    map(x-> addvar!(mp, x), values(var))
#end

#function addvar!(mp::MathPrgm, var::AbstractArray)
#    map(x -> addvar!(mp, x), flatten(var))
#end

#"""
#    addequ!(mp, eqn)
#
#Add an equation to a Mathematical Programm. The argument eqn is a single or list of JuMP.ConstraintRef, not an expression
#"""
#function addequ!(mp::MathPrgm, eqn::JuMP.ConstraintRef, isNL::Bool)
#    push!(mp.equs, (eqn.index.value, isNL))
#end

#function addequ!(mp::MathPrgm, eqn::Vector{ConstraintRef{M,C}}) where {M <: JuMP.AbstractModel, C <: JuMP.AbstractConstraint}
#    map(x -> addequ!(mp, x), eqn)
#end

#function addequ!(mp::MathPrgm, eqn::Vector{JuMP.ConstraintRef})
#    map(x -> addequ!(mp, x), eqn)
#end

#function addequ!(mp::MathPrgm, eqn::AbstractArray, isNL::Bool)
#    map(x -> addequ!(mp, x, isNL), flatten(eqn))
#
#end
#
#function addequ!(mp::MathPrgm, eqn::Vector, isNL::Bool)
#    map(x-> addequ!(mp, x, isNL), values(eqn))
#end


"""
     addovf!(emp, v, args, name, params)

Tag a variable as a supported OVF

# Arguments
- `emp::EMPmaster`: the EMP object
- `v::JuMP.VariableRef`: the variable
- `args::Vector{JuMP.VariableRef}`: the arguments for the OVF problem
- `name::String`: the name of the OVF
- `params::Dict`: the parameters for defining this OVF
"""
function addovf!(emp::EMPmaster, v::JuMP.VariableRef, args::Vector{JuMP.VariableRef}, name::String, params::Dict)
    push!(emp.ovfs, OVF(v, args, name, params))
end

"""
    solveEMP(emp)

Solve an EMP model
"""
function solveEMP(emp::EMPmaster)
    m = getReSHOPModel(emp)

    if (m.mdl == C_NULL)
        m.mdl = ReSHOP.reshop_alloc(m.ctx)
    end

    ReSHOP.emp_init(m.mdl)

    # This is needed to get the metadata initialized
    if (length(emp.mps) > 0)
        ReSHOP.reshop_set_modeltype(m.ctx, ReSHOP.emp)
    end

    # This is needed to get the NL part loaded, before emptree is called
    if emp.backend.nlp_data !== nothing
        MOI.set(emp.backend, MOI.NLPBlock(), JuMP._create_nlp_block_data(emp.backend))
        empty!(emp.backend.nlp_data.nlconstr_duals)
    end

    emp.NLoffset = m.start_nl_cons

    # Deal with the emptree
    emptree(emp, m.mdl)

    # Deal with the OVF functions
    for ovf in emp.ovfs
        ReSHOP.reshop_ovf(m.mdl, ovf)
    end

    optimize!(emp.backend)
end

#########################################################################################################
# End of public API
#########################################################################################################

include("emptree.jl")

function _solveEMP(emp::EMPmaster)

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
        make_equ_index!(emp.equs, emp.backend.internalModel)
    else
        ReSHOP.make_con_index!(reshop_model)
    end

    ctx = ReSHOP.create_reshop_ctx(reshop_model)
    reshop_model.reshop_ctx = ctx

    reshop_emp = ReSHOP.emp_create(ctx)

    if (length(emp.mps) > 0)
        ReSHOP.reshop_set_modeltype(ctx, ReSHOP.emp)
    else
        ReSHOP.reshop_set_modeltype(reshop_model)
    end

    # Define the EMP tree
    emptree(emp, reshop_emp, ctx)

    for ovf in emp.ovfs
        ReSHOP.reshop_ovf(reshop_emp, ovf)
    end

    setvarnames(ctx, emp.backend)

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
        _setvalues(emp.backend, getsolution(emp))
    else
        println("DEBUG: solver failed with status $(reshop_model.solve_exitcode)")
        reshop_model.status = :Error
        reshop_model.solution = fill(NaN,reshop_model.nvar)
        reshop_model.solve_result = "failure"
        reshop_model.solve_result_num = 999
    end

    #    TODO(xhub) move that elsewhere
    # This is left here for now, otherwise, we can't remove the test folder
    # This also plugs a memleak on a context
    ReSHOP.emp_delete(reshop_emp)
end

###################################################################################################
#
#   Accessor and such
#
###################################################################################################

getsolution(emp::EMPmaster) = copy(getReSHOPModel(emp).solution)
getsolvetime(emp::EMPmaster) = getReSHOPModel(emp).solve_time
status(emp::EMPmaster) = status(getReSHOPModel(emp))

function getobjval(mp::MathPrgm)
    reshop_model = getReSHOPModel(mp.emp)
    objvar = ReSHOP.emp_mp_getobjvar(mp.solverobj)
    ReSHOP.ctx_getvarval(reshop_model.reshop_ctx, objvar)
end

# Access to solve results
termination_status(emp::EMPmaster) = JuMP.termination_status(emp.backend)
primal_status(emp::EMPmaster) = JuMP.primal_status(emp.backend)
dual_status(emp::EMPmaster) = JuMP.dual_status(emp.backend)

function Base.show(io::IO, m::EMPmaster)
    println(io, "EMP data structure")
end

end
