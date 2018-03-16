function check_jamsd(m)
    error("Solver $m is not of type JAMSD")
end

# trickery here
# if we store emp directly, print blows up ...
function setemp_solver!(m::JAMSDWriter.JAMSDSolver, emp)
    m.emp = function () return _solveEMP(emp) end
    return true
end

function setemp!(m, emp)
    error("Model $(typeof(m)) not supported: only JuMP and Convex are")
end

function setemp!(m::JuMP.Model, emp)
    setemp_solver!(m.solver, emp)
end

function getJAMSDModel(emp::EMP.Model{JuMP.Model})
    jmodel = emp.model_ds
    if !jmodel.internalModelLoaded
        JuMP.build(emp.model_ds)
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

macro constraintFromExprMP(mp, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @constraint $(esc(mp)).emp.model_ds $dummyconstr )
    quote
        $code
        local gdix = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end

macro constraintFromExprMP2(mp, constr)
    quote
        local cref = @constraint $(esc(mp)).emp.model_ds $(esc(constr))
        local gdix = addequ!($(esc(mp)), cref)
        (cref, gdix)
    end
end

macro NLconstraintFromExprMP(mp, expr)
    dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    code = :( cref = @NLconstraint $(esc(mp)).emp.model_ds $dummyconstr )
    quote
        $code
        local gidx = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
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


