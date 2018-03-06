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

export @variableMP, @objectiveMP, @NLobjectiveMP, @constraintMP, @NLconstraintMP, vipair, NLvipair, solveEMP, MathPrgm, addvar!, addovf!, getsolution, status, get_solve_result, get_solve_result_num, get_model_result, get_model_result_num, get_solve_message, get_solve_exitcode


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


type Model{T<:modeltypes}
    m::T
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    equs::Vector{Tuple{Symbol, Int}}
    ovfs::Vector{OVF}
end

function check_jamsd(m)
    error("Solver $m is not of type JAMSD")
end

# triggery here
# if we store emp directly, print blows up ...
function setemp_solver!(m::JAMSDWriter.JAMSDSolver, emp)
    m.emp = function () return solveEMP(emp) end
    return true
end

function setemp!(m, emp)
    error("Model $(typeof(m)) not supported: only JuMP and Convex are")
end

function setemp!(m::JuMP.Model, emp)
    setemp_solver!(m.solver, emp)
end

function Model(m)
    emp = Model(m, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}(), Vector{Tuple{Symbol, Int}}(), Vector{OVF}())
    setemp!(m, emp)
    return emp
end

function MathPrgm(m::EMP.Model)
    mp = MathPrgm(m, Vector{Int}(), Vector{Int}(), Dict{Int,Int}(), -1, -1, :undef, Vector{MathPrgm}(), Vector{Vector{MathPrgm}}())
    push!(m.mps, mp)
    return mp
end

function Equilibrium(emp::EMP.Model, mps::Vector{MathPrgm})
    push!(emp.equils, mps)
end

function getJAMSDModel(emp::EMP.Model{JuMP.Model})
    jmodel = emp.m
    if !jmodel.internalModelLoaded
        JuMP.build(emp.m)
    end

    return jmodel.internalModel.inner
end

if HAVE_CONVEX
    function getJAMSDModel(emp::EMP.Model{Convex.Problem})
        error("EMP.jl does not yet Convex.jl")
    end
end

function reg_equ(emp::EMP.Model{JuMP.Model}, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.NonlinearExprData}})
    push!(emp.equs, (:NL, cref.idx))
    return length(emp.equs)
end

function reg_equ(emp::EMP.Model{JuMP.Model}, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}}})
    push!(emp.equs, (:quad, cref.idx))
    return length(emp.equs)
end

function reg_equ(emp::EMP.Model{JuMP.Model}, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}})
    push!(emp.equs, (:lin, cref.idx))
    return length(emp.equs)
end

function addvar!(mp::MathPrgm, v::JuMP.Variable)
    push!(mp.vars, v.col)
end

function addvar!(mp::MathPrgm, v::Vector{JuMP.Variable})
    map(x -> addvar!(mp, x), v)
end

function addvar!(mp::MathPrgm, v::JuMP.JuMPArray{JuMP.Variable, 1, Tuple{Int64}})
    map(x-> addvar!(mp, x, m), values(v))
end

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

    local len = length(args)
    if len == 0
        code = :( v = @eval @variable $(mp).emp.m )
    elseif len == 1
        code = :( v = @variable $(mp).emp.m $(esc(args[1])) )
    elseif len == 2
        code = :( v = @eval @variable $(mp).emp.m $(esc(args[2])) $(esc(args[3])) )
    elseif len == 3
        code = :( v = @eval @variable $(mp).emp.m $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) )
    elseif len == 4
        code = :( v = @eval @variable $(mp).emp.m $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) $(esc(args[5])) )
    elseif len == 6
        code = :( v = @eval @variable $(mp).emp.m $(esc(args[2])) $(esc(args[3])) $(esc(args[4])) $(esc(args[5])) $(esc(args[6])) )
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
    code = :( cref = @NLconstraint $(esc(mp)).emp.m $dummyconstr )
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
    code = :( cref = @constraint $(esc(mp)).emp.m $dummyconstr )
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
        cref = @constraint($(esc(mp)).emp.m, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro constraintMP(mp, name, expr)
    quote
        cref = @constraint($(esc(mp)).emp.m, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro NLconstraintMP(mp, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.m, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro NLconstraintMP(mp, name, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.m, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

macro constraintFromExprMP(mp, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @constraint $(esc(mp)).emp.m $dummyconstr )
    quote
        $code
        local gdix = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end

macro constraintFromExprMP2(mp, constr)
    quote
        local cref = @constraint $(esc(mp)).emp.m $(esc(constr))
        local gdix = addequ!($(esc(mp)), cref)
        (cref, gdix)
    end
end

macro NLconstraintFromExprMP(mp, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @NLconstraint $(esc(mp)).emp.m $dummyconstr )
    quote
        $code
        local gidx = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end

function vipair(mp::MathPrgm, expr, var::JuMP.Variable)
#    cref, eidx = @constraintFromExprMP(mp, expr)
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


function addovf!(emp::Model{JuMP.Model}, v::JuMP.Variable, args::Vector{JuMP.Variable}, name::String, params::Dict)
    push!(emp.ovfs, OVF(v, args, name, params))
end



function make_equ_index!(equs, m) error("make_index: unsupported type $(typeof(m))") end

function make_equ_index!(equs::Vector{Tuple{Symbol,Int}}, m::JAMSDWriter.JAMSDNonlinearModel)
    jamsd_model = m.inner

    for elt in Iterators.filter(x -> x[2][1] == :lin, enumerate(equs))
        jamsd_model.nonquad_idx[elt[2][2]] = elt[1]
    end

    start_quad = length(jamsd_model.nonquad_idx)

    for elt in Iterators.filter(x -> x[2][1] == :quad, enumerate(equs))
        jamsd_model.nonquad_idx[elt[2][2]+start_quad] = elt[1]
    end

    start_nl = length(jamsd_model.nonquad_idx)

    for elt in Iterators.filter(x -> x[2][1] == :NL, enumerate(equs))
        jamsd_model.nonquad_idx[elt[2][2]+start_nl] = elt[1]
    end
end

function make_equ_index!(equs::Vector{Tuple{Symbol,Int}}, m::JAMSDWriter.JAMSDLinearQuadraticModel)
    jamsd_model = m.inner

    for elt in Iterators.filter(x -> x[2][1] == :quad, enumerate(equs))
        jamsd_model.quad_idx[elt[2][2]] = elt[1]
    end

    for elt in Iterators.filter(x -> x[2][1] == :lin, enumerate(equs))
        jamsd_model.nonquad_idx[elt[2][2]] = elt[1]
    end

end

function _setvalues(m, x::Vector) error("_setvalues: unsupported type $(typeof(m))") end

function _setvalues(m::JuMP.Model, x::Vector)
    @assert length(m.colVal) == length(x)
    m.colVal[:] = x[:]
end

function setvarnames(ctx, m)
end

function setvarnames(ctx, m::JuMP.Model)
    if !m.internalModel.inner.d.hasvalue
        JAMSDWriter.ctx_setvarnames(ctx, m.colNames)
    end
end

function solveEMP(emp::EMP.Model)
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
        make_equ_index!(emp.equs, emp.m.internalModel)
    else
        JAMSDWriter.make_con_index!(jamsd_model)
    end

    ctx = create_jamsd_ctx(jamsd_model)
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

    setvarnames(ctx, emp.m)
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
        _setvalues(emp.m, getsolution(emp))
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
