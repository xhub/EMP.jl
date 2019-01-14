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

macro constraintFromExprMP(mp, expr)
    #dummyconstr = Expr(:call, esc(:(==)), esc(expr), 0)
    dummyconstr = Expr(:call, :(==), esc(expr), 0)
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
    #dummyconstr = Expr(:call, :(==), esc(expr), 0)
    code = :( cref = @NLconstraint $(esc(mp)).emp.model_ds $dummyconstr )
    quote
        $code
        local gidx = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end
# TODO(xhub) add check warning for an already defined objective equation
"""
    @NLobjectiveMP(mp, sense, expr)

Add a nonlinear objective to a mathematical programm. Sense is either `:Min` or `:Max`.
"""
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

"""
    @objectiveMP(mp, sense, expr)

Add a linear objective to a mathematical programm. Sense is either `:Min` or `:Max`.
"""
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

"""
    @constraintMP(mp, expr)

Add a linear constraint to a mathematical programm
"""
macro constraintMP(mp, expr)
    quote
        cref = @constraint($(esc(mp)).emp.model_ds, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

"""
    @constraintMP(mp, name, expr)

Add a linear constraint (with a identifier `name`) to a mathematical programm
"""
macro constraintMP(mp, name, expr)
    quote
        cref = @constraint($(esc(mp)).emp.model_ds, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

"""
    @NLconstraintMP(mp, expr)

Add a nonlinear constraint to a mathematical programm
"""
macro NLconstraintMP(mp, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.model_ds, $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

"""
    @NLconstraintMP(mp, name, expr)

Add a nonlinear constraint (with a identifier `name`) to a mathematical programm
"""
macro NLconstraintMP(mp, name, expr)
    quote
        cref = @NLconstraint($(esc(mp)).emp.model_ds, $(esc(name)), $(esc(expr)))
        gdix = addequ!($(esc(mp)), cref)
        cref
    end
end

"""
    vipair(mp, expr, var)

Add an affine variational inequality relationship between the variable `var` and the mapping `expr`

If there is no additional constraint of the feasible set of `var` besides lower and upper bounds,
then this is defines a Mixed Complementarity Problem
    expr ⟂ lb ≤ var ≤ ub
Or
    0 ∈ expr + N_[lb, ub] (var)
"""
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

macro vipair(mp, expr, var)
    quote
        return vipair($(esc(mp)), $(esc(expr)), $(esc(var)))
    end
end

"""
    @NLvipair(mp, expr, var)

Add a nonlinear variational inequality relationship between the variable `var` and the mapping `expr`

If there is no additional constraint of the feasible set of `var` besides lower and upper bounds,
then this is defines a Mixed Complementarity Problem
    expr ⟂ lb ≤ var ≤ ub
Or
    0 ∈ expr + N_[lb, ub] (var)
"""
macro NLvipair(mp, expr, var::JuMP.Variable)
    cref, eidx = @NLconstraintFromExprMP(mp, expr)
    mp.matching[var.col] = eidx
    cref
end

#function NLvipair(mp::MathPrgm, expr::Vector{Any}, var::Vector{JuMP.Variable})
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        NLvipair(mp, expr[i], var[i])
#    end
#end


