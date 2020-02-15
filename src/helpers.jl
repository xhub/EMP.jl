function check_reshop(m)
    error("Solver $m is not of type ReSHOP")
end

# trickery here
# if we store emp directly, print blows up ...
function setemp_solver!(m::ReSHOP.Optimizer, emp)
#    ReSHOP.reshop_setemp!(m, function () return _solveEMP(emp) end)
    return true
end

function setemp!(m, emp)
    error("Model $(typeof(m)) not supported: only JuMP and Convex are")
end

function setemp!(m::JuMP.Model, emp)
    setemp_solver!(m.moi_backend.optimizer.model, emp)
end

function getReSHOPModel(emp::EMP.Model)
    jmodel = emp.model_ds

    return jmodel.moi_backend.optimizer.model
end

function _setvalues(m, x::Vector) error("_setvalues: unsupported type $(typeof(m))") end

function _setvalues(m::JuMP.Model, x::Vector)
    @assert length(m.colVal) == length(x)
    m.colVal[:] = x[:]
end

function setvarnames(ctx, m)
end

function setvarnames(ctx, m::JuMP.Model)
    if isnothing(m.internalModel.inner.d)
        ReSHOP.ctx_setvarnames(ctx, m.colNames)
    end
end


