__precompile__()

module EMP

using JuMP
using JAMSDWriter

HAVE_CONVEX = try using Convex; true catch false end

if HAVE_CONVEX
    modeltypes = Union{JuMP.Model,Convex.Problem}
else
    modeltypes = JuMP.Model
end

export @variableMP, @objectiveMP, @NLobjectiveMP, @constraintMP, @NLconstraintMP, vipair, NLvipair, solveEMP, _solveEMP, MathPrgm, addvar!, addequ!, addovf!, getsolution, status, get_solve_result, get_solve_result_num, get_model_result, get_model_result_num, get_solve_message, get_solve_exitcode, solve

" Mathematical Programm representation "
type MathPrgm
    emp::Any
    vars::Vector{Int}
    equs::Vector{Int}
    matching::Dict{Int,Int}
    objequ::Int
    objvar::Int
    sense::Symbol
    mp::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
end

" Optimal Value Function (OVF) representation "
type OVF
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
type Model{T<:modeltypes}
    model_ds::T
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    equs::Vector{Tuple{Symbol, Int}}
    ovfs::Vector{OVF}
end

"""
    Model([modeling_pkg = "JuMP", solver = "jams"])

Create an EMP master object, using the modeling package given as argument to store the variables and equations.
The solver argument is used to construct the JAMSDSolver object.
"""
function Model(; modeling_pkg = "JuMP", solver = "jams")
    if modeling_pkg == "JuMP"
        solver_jamsd = JAMSDWriter.JAMSDSolver(solver)
        model_ds = JuMP.Model(solver=solver_jamsd)
        emp = Model(model_ds, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), Vector{Tuple{Symbol, Int}}(), Vector{OVF}())
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
function Model{T<:modeltypes}(model_ds::T)
    emp = Model(model_ds, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), Vector{Tuple{Symbol, Int}}(), Vector{OVF}())

    setemp!(model_ds, emp)

    return emp
end

"""
    MathPrgm(m::EMP.Model)

Create a Mathematical Programm in the EMP master object
"""
function MathPrgm(m::EMP.Model)
    mp = MathPrgm(m, Vector{Int}(), Vector{Int}(), Dict{Int,Int}(), -1, -1, :undef, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}())
    push!(m.mps, mp)
    return mp
end

" Create an Equilibrium object (not implemented) "
function Equilibrium(emp::EMP.Model, mps::Vector{MathPrgm})
    push!(emp.equils, mps)
end

include("helpers.jl")

"""
    addvar!(mp, var)

Add a variable to a Mathematical Programm
"""
function addvar!(mp::MathPrgm, v::JuMP.Variable)
    push!(mp.vars, v.col)
end

function addvar!(mp::MathPrgm, v::Vector{JuMP.Variable})
    map(x -> addvar!(mp, x), v)
end

function addvar!(mp::MathPrgm, v::JuMP.JuMPArray{JuMP.Variable, 1, Tuple{Int64}})
    map(x-> addvar!(mp, x, m), values(v))
end

"""
    addequ!(mp, eqn)

Add an equation to a Mathematical Programm
"""
function addequ!(mp::MathPrgm, cref::ConstraintRef)
    gidx = reg_equ(mp.emp, cref)
    push!(mp.equs, gidx)
    return gidx
end

function addequ!{M<:JuMP.AbstractModel,C<:JuMP.AbstractConstraint}(mp::MathPrgm, cref::Vector{ConstraintRef{M,C}})
    map(x -> addequ!(mp, x), cref)
end

function addequ!(mp::MathPrgm, cref::Vector{JuMP.ConstraintRef})
    map(x -> addequ!(mp, x), cref)
end

function addequ!(mp::MathPrgm, cref::JuMP.JuMPArray{JuMP.ConstraintRef, 1, Tuple{Int64}})
    map(x-> addequ!(mp, x), values(cref))
end

"""
Add a variable to a Mathematical Programm. See JuMP `@variable` for examples
"""
macro variableMP(mp, args...)
    mp = esc(mp)
    # Pick out keyword arguments
    if Meta.isexpr(args[1], :parameters) # these come if using a semicolon
        kwargs = args[1]
        args = args[2:end]
    else
        kwargs = Expr(:parameters)
    end
    kwsymbol = VERSION < v"0.6.0-dev.1934" ? :kw : :(=) # changed by julia PR #19868
    append!(kwargs.args, collect(Iterators.filter(x -> Meta.isexpr(x, kwsymbol), collect(args)))) # comma separated
    args = collect(Iterators.filter(x -> !Meta.isexpr(x, kwsymbol), collect(args)))
    if length(kwargs.args) > 0
        error("Please do not use keyword argument in @variableMP")
    end

    # TODO(xhub) fix this with proper syntax
    # dummyconstr = Expr(:call, @variable, $(mp).emp.m, $(esc(args[1])))
    #then extend the arguments
    local len = length(args)
    if len == 0
        code = :( v = @eval @variable $(mp).emp.model_ds )
    elseif len == 1
        code = :( v = @variable $(mp).emp.model_ds $(esc(args[1])) )
    elseif len == 2
        code = :( v = @eval @variable $(mp).emp.model_ds $(esc(args[2])) $(esc(args[3])) )
    elseif len == 3
        code = :( v = @eval @variable $(mp).emp.model_ds $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) )
    elseif len == 4
        code = :( v = @eval @variable $(mp).emp.model_ds $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) $(esc(args[5])) )
    elseif len == 6
        code = :( v = @eval @variable $(mp).emp.model_ds $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) $(esc(args[5])) $(esc(args[6])) )
    end
    quote
        if $len > 6
            error("unsupported syntax $(esc(args[2]))")
        end
        $code
        addvar!($(mp), v)
        v
    end
end

macro NLobjectiveMP(mp, sense, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @NLconstraint $(esc(mp)).emp.model_ds $dummyconstr )
    quote
        $code
        gidx = addequ!($(esc(mp)), cref)
        $(esc(mp)).objequ = gidx
        $(esc(mp)).sense = $(esc(sense))
        cref
    end
end

macro objectiveMP(mp, sense, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @constraint $(esc(mp)).emp.model_ds $dummyconstr )
    quote
        $code
        gidx = addequ!($(esc(mp)), cref)
        $(esc(mp)).objequ = gidx
        $(esc(mp)).sense = $(esc(sense))
        cref
    end
end

macro constraintMP(mp, expr)
    quote
        cref = @constraint($(esc(mp)).emp.model_ds, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro constraintMP(mp, name, expr)
    quote
        cref = @constraint($(esc(mp)).emp.model_ds, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro NLconstraintMP(mp, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.model_ds, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro NLconstraintMP(mp, name, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.model_ds, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

"""
    vipair(mp, expr, var)

Add a variational inequality relationship between the variable `var` and the mapping `expr`

If there is no additional constraint of the feasible set of `var` besides lower and upper bouds,
then this is defines a Mixed Complementarity Problem
    expr ⟂ lb ≤ var ≤ ub
Or
    0 ∈ expr + N_[lb, ub] (var)
"""
function vipair(mp::MathPrgm, expr, var::JuMP.Variable)
    cref, eidx = @constraintFromExprMP2(mp, expr == 0.)
    mp.matching[var.col] = eidx
    cref
end

function vipair(mp::MathPrgm, expr::Vector, var::Vector{JuMP.Variable})
    @assert length(expr) == length(var)
    for i in 1:length(var)
        vipair(mp, expr[i], var[i])
    end
end

macro NLvipair(mp, expr, var)
    cref, eidx = @NLconstraintFromExprMP(mp, expr)
    mp.matching[var.col] = eidx
    cref
end

function NLvipair(mp::MathPrgm, expr::Vector{Any}, var::Vector{JuMP.Variable})
    @assert length(expr) == length(var)
    for i in 1:length(var)
        NLvipair(mp, expr[i], var[i])
    end
end

"""
Tag a variable as a supported OVF

# Arguments
- emp: the EMP object
- v: the variable
- args: the arguments for the OVF problem
- `name::String`: the name of the OVF
- `params::Dict`: the parameters for defining this OVF

"""
function addovf!(emp::Model{JuMP.Model}, v::JuMP.Variable, args::Vector{JuMP.Variable}, name::String, params::Dict)
    push!(emp.ovfs, OVF(v, args, name, params))
end

"""
    solveEMP(emp)

Solve an EMP model
"""
function solveEMP{T<:modeltypes}(emp::EMP.Model{T})
    error("solve not implemented for a model stat structure of type $(typeof(emp.model_ds))")
end

function solveEMP(emp::EMP.Model{JuMP.Model})
    JuMP.solve(emp.model_ds)
end

#########################################################################################################
# End of public API
#########################################################################################################

function _solveEMP(emp::EMP.Model)
    # Steps
    # 1. Define the full models and agents
    # 2. Solve
    # 3. Report values
    jamsd_model = getJAMSDModel(emp)

    JAMSDWriter.make_var_index!(jamsd_model)
#    JAMSDWriter.make_con_index!(jamsd_model)

    ########
    #
    # We have to mirror the behavior of JAMSDWriter, where the nonlinear equation are put before
    # the linear ones.
    #
    # TODO(xhub) fix this mess ...
    #
    #####

    if (length(emp.mps) > 0)
        make_equ_index!(emp.equs, emp.model_ds.internalModel)
    else
        JAMSDWriter.make_con_index!(jamsd_model)
    end

    ctx = JAMSDWriter.create_jamsd_ctx(jamsd_model)
    jamsd_model.jamsd_ctx = ctx

    jamsd_emp = JAMSDWriter.emp_create(ctx)

    JAMSDWriter.emp_mp_ensure(jamsd_emp, length(emp.mps))

    for (idx, mp) in enumerate(emp.mps)
        jamsd_mp = JAMSDWriter.jamsd_declare_mathprgm(mp, ctx, idx)
        JAMSDWriter.emp_mp_store(jamsd_emp, jamsd_mp)
    end

    for ovf in emp.ovfs
        JAMSDWriter.jamsd_ovf(jamsd_emp, ovf)
    end

    if (length(emp.mps) > 0)
        JAMSDWriter.jamsd_set_modeltype(ctx, JAMSDWriter.emp)
    else
        JAMSDWriter.jamsd_set_modeltype(jamsd_model)
    end

    setvarnames(ctx, emp.model_ds)
#    JAMSDWriter.emp_mp_to_agent(ctx, jamsd_emp)

    gams_ctx = JAMSDWriter.jamsd_setup_gams()
    jamsd_model.jamsd_ctx_dest = gams_ctx

    jamsd_model.solve_exitcode = JAMSDWriter.jamsd_solve(ctx, gams_ctx, jamsd_model.solver_name, jamsd_emp)

    ###############################################################################################
    #
    #  Report values and codes
    #
    ###############################################################################################

    if jamsd_model.solve_exitcode == 0
        JAMSDWriter.emp_report_values(jamsd_emp, ctx)
        JAMSDWriter.report_results_common(jamsd_model)
        _setvalues(emp.model_ds, getsolution(emp))
    else
        println("DEBUG: solver failed with status $(jamsd_model.solve_exitcode)")
        jamsd_model.status = :Error
        jamsd_model.solution = fill(NaN,jamsd_model.nvar)
        jamsd_model.solve_result = "failure"
        jamsd_model.solve_result_num = 999
    end

    JAMSDWriter.emp_delete(jamsd_emp)
end

###################################################################################################
#
#   Accessor and such
#
###################################################################################################

getsolution(emp::EMP.Model) = copy(getJAMSDModel(emp).solution)
getsolvetime(emp::EMP.Model) = getJAMSDModel(emp).solve_time
status(emp::EMP.Model) = status(getJAMSDModel(emp))

# Access to solve results
get_solve_result(emp::EMP.Model) = getJAMSDModel(emp).solve_result
get_solve_result_num(emp::EMP.Model) = getJAMSDModel(emp).solve_result_num
get_model_result(emp::EMP.Model) = getJAMSDModel(emp).model_result
get_model_result_num(emp::EMP.Model) = getJAMSDModel(emp).model_result_num
get_solve_message(emp::EMP.Model) = getJAMSDModel(emp).solve_message
get_solve_exitcode(emp::EMP.Model) = getJAMSDModel(emp).solve_exitcode

end
