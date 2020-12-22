
# Simple example of JuMP extension used in the tests to check that JuMP works well with extensions
# The main difference between `JuMP.Model` and `JuMPExtension.MathPrgm` is the fact that in `add_variable` (resp. `add_constraint`),
# `JuMP.Model` applies the modification to its `moi_backend` field while
# `JuMPExtension.MathPrgm` stores the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

using JuMP

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


struct ConstraintIndex
    value::Int # Index in `model.constraints`
end

" Mathematical Programm representation "
mutable struct MathPrgm <: JuMP.AbstractModel
    emp_master
    vars::Vector{Int}
    equs::Vector{Tuple{Int, Bool}}
    matching::Dict{Int,Tuple{Int,Bool}}
    objectivesense::Union{Symbol,MOI.OptimizationSense}
    mps::Vector{MathPrgm}
    equils::Vector{Vector{MathPrgm}}
    solverobj::Ptr{}

    objective_function::JuMP.AbstractJuMPScalar
    obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
end


Base.broadcastable(model::MathPrgm) = Ref(model)

JuMP.object_dictionary(model::MathPrgm) = model.obj_dict

# Variables
#Base.copy(v::MyVariableRef) = v
#Base.copy(v::MyVariableRef, new_model::MathPrgm) = MyVariableRef(new_model, v.idx)

#Base.:(==)(v::MyVariableRef, w::MyVariableRef) = v.model === w.model && v.idx == w.idx
#Base.broadcastable(v::MyVariableRef) = Ref(v)
#JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w
JuMP.variable_type(::MathPrgm) = JuMP.VariableRef
function JuMP.add_variable(m::MathPrgm, v::JuMP.AbstractVariable, name::String="")
    vref = JuMP.add_variable(model.emp.backend, v, name)
    push!(mp.vars, vref.index.value)
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
    remove!(model.emp.backend, vref.index.value)
    delete(model.emp.backend, vref)
end
function JuMP.delete(model::MathPrgm, vrefs::Vector{JuMP.VariableRef})
    JuMP.delete.(model, vrefs)
end
function JuMP.is_valid(model::MathPrgm, vref::JuMP.VariableRef)
    return (model === vref.model &&
            vref.idx in keys(model.variables))
end
JuMP.num_variables(m::MathPrgm) = length(m.vars)

# Constraints
JuMP.constraint_type(::MathPrgm) = JuMP.ConstraintRef
function JuMP.add_constraint(model::MathPrgm, c::JuMP.AbstractConstraint,
                             name::String="")
    cref = JuMP.add_constraint(model.emp.backend, c, name)
    push!(mp.equs, cref.index.value)
    cref
end
function JuMP.delete(model::MathPrgm, cref::JuMP.ConstraintRef)
    @assert JuMP.is_valid(model, constraint_ref)
    remove!(model.equs, cref.index.value)
    delete(model.emp.backend, cref)
end
function JuMP.delete(model::MathPrgm, con_refs::Vector{<:JuMP.ConstraintRef})
    JuMP.delete.(model, con_refs)
end
function JuMP.is_valid(model::MathPrgm, constraint_ref::JuMP.ConstraintRef)
    return (model === constraint_ref.model &&
            constraint_ref.index in keys(model.constraints))
end
function JuMP.num_constraints(model::MathPrgm,
    F::Type{<:JuMP.AbstractJuMPScalar},
    S::Type{<:MOI.AbstractSet})
    crefs = all_constraints(model.emp.backend, F, S)
    return count(cref -> cref.index.value in model.equs, crefs)
end
function JuMP.num_constraints(model::MathPrgm,
    T::Type{<:Vector{F}},
    S::Type{<:MOI.AbstractSet}) where F<:JuMP.AbstractJuMPScalar
    crefs = all_constraints(model.emp.backend, T, S)
    return count(cref -> cref.index.value in model.equs, crefs)
end


# Objective
function JuMP.set_objective_function(model::MathPrgm, f::JuMP.AbstractJuMPScalar)
    m.objective_function = f
end
function JuMP.set_objective_function(m::MathPrgm, f::Real)
    m.objective_function = JuMP.GenericAffExpr{Float64, JuMP.VariableRef}(f)
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
    n = length(model.constraints)
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
