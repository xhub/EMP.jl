
# Simple example of JuMP extension used in the tests to check that JuMP works well with extensions
# The main difference between `JuMP.Model` and `JuMPExtension.MathPrgm` is the fact that in `add_variable` (resp. `add_constraint`),
# `JuMP.Model` applies the modification to its `moi_backend` field while
# `JuMPExtension.MathPrgm` stores the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

using JuMP

function remove!(a, item)
    deleteat!(a, findall(x -> x==item, a))
end


struct ConstraintIndex
    value::Int # Index in `model.constraints`
end

struct FakeNLP
  nlconstr::Vector{Any}
  indices::AbstractVector{Int}
end

" Mathematical Programm representation "
mutable struct MathPrgm <: JuMP.AbstractModel
    emp
    vars::Vector{Int}
    cons::Vector{MOI.ConstraintIndex} # constraint references
    equs::Vector{MOI.ConstraintIndex} # equation (most likely VI mapping) references
    objectivesense::Union{Symbol,MOI.OptimizationSense}
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    solverobj::Ptr{}

    objective_function::Union{MOI.AbstractScalarFunction, Int}
    obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
    nlp_data
    function MathPrgm(m)
        mp = new(m,
                 Vector{Int}(),                # vars
                 Vector{MOI.ConstraintIndex}(),   # cons
                 Vector{MOI.ConstraintIndex}(),   # equs
                 MOI.FEASIBILITY_SENSE,        # objectivesense
                 Vector{MathPrgm}(),           # mp
                 Vector{Vector{MathPrgm}}(),   # equils
                 C_NULL,
                 zero(MOI.ScalarQuadraticFunction{Float64}),
                 Dict{Symbol, Any}(),
                 nothing)
        push!(m.mps, mp)
        return mp
    end
end


Base.broadcastable(model::MathPrgm) = Ref(model)

JuMP.object_dictionary(model::MathPrgm) = model.obj_dict

# Helpers

isNL(c::JuMP.AbstractConstraint) = false
isNL(c::JuMP._NonlinearConstraint) = true

function _add_ei(model::MathPrgm, cref::JuMP.ConstraintRef, set::MOI.AbstractSet)
  push!(model.cons, cref.index)
end

function _add_ei(model::MathPrgm, cref::JuMP.ConstraintRef, set::MOI.Complements)
  push!(model.equs, cref.index)
end

function add_ei(model::MathPrgm, cref::JuMP.ConstraintRef, c::Union{JuMP.VectorConstraint, JuMP.ScalarConstraint})
  _add_ei(model, cref, c.set)
end

function add_ei(model::MathPrgm, cref::JuMP.ConstraintRef, c::JuMP.AbstractConstraint)
  push!(model.cons, cref.index)
end


# Variables
#Base.copy(v::MyVariableRef) = v
#Base.copy(v::MyVariableRef, new_model::MathPrgm) = MyVariableRef(new_model, v.idx)

#Base.:(==)(v::MyVariableRef, w::MyVariableRef) = v.model === w.model && v.idx == w.idx
#Base.broadcastable(v::MyVariableRef) = Ref(v)
#JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w
JuMP.variable_type(::MathPrgm) = JuMP.VariableRef

# This function does the real add
function JuMP.add_variable(model::MathPrgm, v::JuMP.AbstractVariable, name::String="")
    vref = JuMP.add_variable(model.emp.backend, v, name)
    push!(model.vars, vref.index.value)
    vref
end

function JuMP.add_variable(model::MathPrgm, variable::JuMP.VariableConstrainedOnCreation, name::String)
    var_ref = JuMP.add_variable(model, variable.scalar_variable, name)
    JuMP.add_constraint(model, JuMP.ScalarConstraint(var_ref, variable.set))
    return var_ref
end

function JuMP.add_variable(model::MathPrgm, variable::JuMP.VariablesConstrainedOnCreation, names)
    var_refs = JuMP.add_variable.(model, variable.scalar_variables,
                                  JuMP.vectorize(names, variable.shape))
    JuMP.add_constraint(model, JuMP.VectorConstraint(var_refs, variable.set))
    return JuMP.reshape_vector(var_refs, variable.shape)
end

function JuMP.delete(model::MathPrgm, vref::JuMP.VariableRef)
    @assert JuMP.is_valid(model.emp.backend, vref)
    remove!(model.vars, vref.index.value)
    delete(model.emp.backend, vref)
end

function JuMP.delete(model::MathPrgm, vrefs::Vector{JuMP.VariableRef})
    JuMP.delete.(model, vrefs)
end
JuMP.is_valid(model::MathPrgm, vref::JuMP.VariableRef) = JuMP.is_valid(model.emp.backend, vref)
JuMP.num_variables(model::MathPrgm) = length(model.vars)

# Constraints
JuMP.constraint_type(::MathPrgm) = JuMP.ConstraintRef
function JuMP.add_constraint(model::MathPrgm, c::JuMP.AbstractConstraint,
                             name::String="")
    cref = JuMP.add_constraint(model.emp.backend, c, name)
    add_ei(model, cref, c)
    cref
end

function JuMP.delete(model::MathPrgm, cref::JuMP.ConstraintRef)
    @assert JuMP.is_valid(model, cref)
    remove!(model.equs, cref.index)
    delete(model.emp.backend, cref)
end
function JuMP.delete(model::MathPrgm, crefs::Vector{<:JuMP.ConstraintRef})
    JuMP.delete.(model, crefs)
end
function JuMP.is_valid(model::MathPrgm, cref::JuMP.ConstraintRef)
    return (model === cref.model &&
            (cref.index in keys(model.cons) || cref.index in keys(model.equs)))
end
function JuMP.num_constraints(model::MathPrgm,
    F::Type{<:JuMP.AbstractJuMPScalar},
    S::Type{<:MOI.AbstractSet})
    crefs = all_constraints(model.emp.backend, F, S)
    return count(cref -> cref.index in model.equs, crefs) + count(cref -> cref.index in model.cons, crefs)
end
function JuMP.num_constraints(model::MathPrgm,
    T::Type{<:Vector{F}},
    S::Type{<:MOI.AbstractSet}) where F<:JuMP.AbstractJuMPScalar
    crefs = all_constraints(model.emp.backend, T, S)
    return count(cref -> cref.index in model.equs, crefs) + count(cref -> cref.index in model.cons, crefs)
end


# Objective
function JuMP.set_objective_function(model::MathPrgm, f::JuMP.AbstractJuMPScalar)
  model.objective_function = JuMP.moi_function(f)
  return
end

function JuMP.set_objective_function(model::MathPrgm, f::Real)
  model.objective_function = JuMP.moi_function(JuMP.GenericAffExpr{Float64, JuMP.VariableRef}(f))
  return
end

function JuMP.set_objective_function(model::MathPrgm, f::JuMP._NonlinearExprData)
    jump_model = model.emp.backend
    JuMP._init_NLP(jump_model)
    nlcons = JuMP._NonlinearConstraint(f, 0., 0.)
    push!(jump_model.nlp_data.nlconstr, nlcons)
    model.objective_function = length(jump_model.nlp_data.nlconstr)
    return
end

function JuMP._init_NLP(model::MathPrgm)
    jump_model = model.emp.backend
    JuMP._init_NLP(jump_model)
    model.nlp_data = FakeNLP(Vector{Any}(), Int[])
    return
end

JuMP.objective_sense(model::MathPrgm) = model.objectivesense
function JuMP.set_objective_sense(model::MathPrgm, sense)
    model.objectivesense = sense
end
JuMP.objective_function_type(model::MathPrgm) = typeof(model.objective_function)
JuMP.objective_function(model::MathPrgm) = model.objective_function
function JuMP.objective_function(model::MathPrgm, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(model.objective_function isa FT)
        throw(InexactError(:objective_function, FT,
                           typeof(model.objective_function)))
    end
    return model.objective_function::FT
end

# Names
function JuMP.variable_by_name(model::MathPrgm, name::String)
    return variable_by_name(model.emp.backend, name)
end
function JuMP.constraint_by_name(model::MathPrgm, name::String)
    return constraint_by_name(model.emp.backend, name)
end

# Show
function JuMP.show_backend_summary(io::IO, model::MathPrgm) end
function JuMP.show_objective_function_summary(io::IO, model::MathPrgm)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(model))
end
function JuMP.objective_function_string(print_mode, model::MathPrgm)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end
_plural(n) = (isone(n) ? "" : "s")
function JuMP.show_constraints_summary(io::IO, model::MathPrgm)
    n = length(model.equs)
    print(io, "Constraint", _plural(n), ": ", n)
end
function JuMP.constraints_string(print_mode, model::MathPrgm)
    strings = String[]
    # Sort by creation order, i.e. ConstraintIndex value
    crefs = all_constraints(model.emp.backend)
    constraints = sort(crefs, by = c -> c.first.index.value)
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    return strings
end

JuMP.all_variables(model::MathPrgm) = JuMP.all_variables(model.emp.backend)
#########################################################################
# EMP add
#########################################################################

JuMP._parse_NL_expr_runtime(model::MathPrgm, a::Any, b::Any, c::Any, d::Any) = JuMP._parse_NL_expr_runtime(model.emp.backend, a, b, c, d)
