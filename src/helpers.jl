function check_reshop(m)
    error("Solver $m is not of type ReSHOP")
end

# trickery here
# if we store emp directly, print blows up ...
function setemp_solver!(m::ReSHOP.ReSHOPSolver, emp)
    ReSHOP.reshop_setemp!(m, function () return _solveEMP(emp) end)
    return true
end

function setemp!(m, emp)
    error("Model $(typeof(m)) not supported: only JuMP and Convex are")
end

function setemp!(m::JuMP.Model, emp)
    setemp_solver!(m.solver, emp)
end

function getReSHOPModel(emp::EMP.Model)
    jmodel = emp.model_ds
    if !jmodel.internalModelLoaded
        JuMP.build(emp.model_ds)
    end

    return jmodel.internalModel.inner
end

function reg_equ(emp::EMP.Model, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.NonlinearExprData}})
    push!(emp.equs, (:NL, cref.idx))
    return length(emp.equs)
end

function reg_equ(emp::EMP.Model, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}}})
    push!(emp.equs, (:quad, cref.idx))
    return length(emp.equs)
end

function reg_equ(emp::EMP.Model, cref::JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}})
    push!(emp.equs, (:lin, cref.idx))
    return length(emp.equs)
end

function make_equ_index!(equs, m) error("make_index: unsupported type $(typeof(m))") end

function make_equ_index!(equs::Vector{Tuple{Symbol,Int}}, m::ReSHOP.ReSHOPNonlinearModel)
    reshop_model = m.inner

    for elt in Iterators.filter(x -> x[2][1] == :lin, enumerate(equs))
        reshop_model.nonquad_idx[elt[2][2]] = elt[1]
    end

    start_quad = length(reshop_model.nonquad_idx)

    for elt in Iterators.filter(x -> x[2][1] == :quad, enumerate(equs))
        reshop_model.nonquad_idx[elt[2][2]+start_quad] = elt[1]
    end

    start_nl = length(reshop_model.nonquad_idx)

    for elt in Iterators.filter(x -> x[2][1] == :NL, enumerate(equs))
        reshop_model.nonquad_idx[elt[2][2]+start_nl] = elt[1]
    end
end

function make_equ_index!(equs::Vector{Tuple{Symbol,Int}}, m::ReSHOP.ReSHOPLinearQuadraticModel)
    reshop_model = m.inner

    for elt in Iterators.filter(x -> x[2][1] == :quad, enumerate(equs))
        reshop_model.quad_idx[elt[2][2]] = elt[1]
    end

    for elt in Iterators.filter(x -> x[2][1] == :lin, enumerate(equs))
        reshop_model.nonquad_idx[elt[2][2]] = elt[1]
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
        ReSHOP.ctx_setvarnames(ctx, m.colNames)
    end
end


