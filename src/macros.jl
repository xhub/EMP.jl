import Base.Iterators: flatten

function ast_sym_esc!(array_args)
    for i=1:length(array_args)
        if (typeof(array_args[i]) == Symbol)
            array_args[i] = esc(array_args[i])
        elseif (typeof(array_args[i]) == Expr)
            ast_sym_esc!(array_args[i].args)
        end
    end
end

function ast_sym2_esc!(ast_node, skipsym=Vector{Symbol}(undef,0))
#    println("ast node is $(ast_node)")
    if (ast_node.head == :call || ast_node.head == :callmacro)
        s=2
    elseif (ast_node.head == :generator)
        # we have faith here that the iterator has a structure like
        # ast_node.args[2] ~ i=1:/N
#        println("skipping generator var $(ast_node.args[2].args[1])")
        push!(skipsym, ast_node.args[2].args[1])
#        println("new skipsym $(skipsym)")
        # now visit ast_node.args[1] and ast_node.args[2].args[2] (i=1:N)
        ast_sym2_esc!(ast_node.args[1], skipsym)
        ast_sym2_esc!(ast_node.args[2].args[2], skipsym)
        return
    else
        s=1
    end
    for i=s:length(ast_node.args)
#        println("skipsym is $(skipsym[:])")
        if (typeof(ast_node.args[i]) == Symbol && !(ast_node.args[i] in skipsym))
#            println("escaping $(ast_node.args[i])")
            ast_node.args[i] = esc(ast_node.args[i])
        elseif (typeof(ast_node.args[i]) == Expr)
            ast_sym2_esc!(ast_node.args[i], skipsym)
        end
    end
end

function ast_args_esc!(array_args)
    for i=1:length(array_args)
        if (typeof(array_args[i]) == Expr)
            array_args[i] = esc(array_args[i])
        end
    end
end

"""
Add a variable to a Mathematical Programm. See JuMP `@variable` for examples
"""
macro variableMP(mp, args...)
#    mp = esc(mp)
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
    len = length(args)

    #ar1 = :($(esc(args[1])))
    dump(args)

    ar1 = args[1]
    if (typeof(ar1) == Symbol)
        varname = esc(ar1)
    elseif (ar1.head == :ref)
        # This is the case there we have var[...]

        varname = esc(ar1.args[1])
        ast_sum_esc!(ar1.args)
    elseif (ar1.head == :call || ar1.head == :comparison)
        # This is the case there we have var >= 0 or var[...] >= 0
        if (ar1.head == :comparison)
            posvar = 3
        else
            if (length(ar1.args) == 3)
                posvar = 2
            else
                error("Could not parse variable declaration $(ar1); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
            end
        end

        if (typeof(ar1.args[posvar]) == Symbol)
            varname = esc(ar1.args[posvar])
        elseif (typeof(ar1.args[posvar]) == Expr)
            if (typeof(ar1.args[posvar].args[1]) != Symbol)
                error("Could not parse variable declaration $(ar1); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
            end
            varname = esc(ar1.args[posvar].args[1])

            #now escape some variables in the indices
            ast_sym_esc!(ar1.args[posvar].args[2:end])
        else
            error("Could not parse variable declaration $(ar1); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
        end
    else
        error("Could not parse variable declaration $(ar1); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
    end

    dump(args)
    return quote
        mmp = $(esc(mp)).emp.model_ds
        #varname = $(esc(args[1]))
        if $len > 6
            error("unsupported syntax $(esc(args[2]))")
        end
        # TODO: use eval(Expr(:call, variable, ...))
        $varname = @variable mmp $args
        addvar!($(esc(mp)), $varname)
        $varname
    end
end

macro constraintFromExprMP(mp, expr)
    return quote
        cref = @constraint $(esc(mp)).emp.model_ds $expr == 0
        gdix = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end

macro NLconstraintFromExprMP(mp, expr)
    return quote
        cref = @NLconstraint $(esc(mp)).emp.model_ds $expr == 0
        gidx = addequ!($(esc(mp)), cref)
        (cref, gidx)
    end
end

# TODO(xhub) add check warning for an already defined objective equation
"""
    @NLobjectiveMP(mp, sense, expr)

Add a nonlinear objective to a mathematical programm. Sense is either `:Min` or `:Max`.
"""
macro NLobjectiveMP(mp, sense, expr)
    ast_sym2_esc!(expr)
    dummyconstr = Expr(:call, :(==), expr, 0)
    mmp = esc(mp)
    return quote
        model = $(mmp).emp.model_ds
        cref = @NLconstraint model $dummyconstr
        gidx = addequ!($mmp, cref)
        $(mmp).objequ = gidx
        $(mmp).sense = $(esc(sense))
        cref
    end
end

"""
    @objectiveMP(mp, sense, expr)

Add a linear objective to a mathematical programm. Sense is either `:Min` or `:Max`.
"""
macro objectiveMP(mp, sense, expr)
    return quote
        mmp = $(esc(mp)).emp.model_ds
        cref = @constraint mmp $expr == 0
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
    ast_args_esc!(expr.args)
    return quote
        mmp = $(esc(mp))
        model = mmp.emp.model_ds
        cref = @constraint model $expr
        gdix = addequ!(mmp, cref)
        cref
    end
end

"""
    @constraintMP(mp, name, expr)

Add a linear constraint (with a identifier `name`) to a mathematical programm
"""
macro constraintMP(mp, name, expr)
    ast_args_esc!(expr.args)
    return quote
        mmp = $(esc(mp))
        model = mmp.emp.model_ds
        cref = @constraint model $name $expr
        gdix = addequ!(mmp, cref)
        cref
    end
end

"""
    @NLconstraintMP(mp, expr)

Add a nonlinear constraint to a mathematical programm
"""
macro NLconstraintMP(mp, expr)
    ast_sym2_esc!(expr)
    return quote
        mmp = $(esc(mp))
        model = mmp.emp.model_ds
        cref = @NLconstraint model $expr
        gdix = addequ!(mmp, cref)
        cref
    end
end

"""
    @NLconstraintMP(mp, name, expr)

Add a nonlinear constraint (with a identifier `name`) to a mathematical programm
"""
macro NLconstraintMP(mp, name, expr)
    ast_sym2_esc!(expr)
    return quote
        mmp = $(esc(mp))
        model = mmp.emp.model_ds
        cref = @NLconstraint model $name $expr
        gdix = addequ!(mmp, cref)
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
#function vipair(mp::MathPrgm, expr, var::JuMP.Variable)
#    cref, eidx = @constraintFromExprMP(mp, expr)
#    mp.matching[var.col] = eidx
#    cref
#end
#
#function vipair(mp::MathPrgm, expr::Vector, var::Vector{JuMP.Variable})
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        vipair(mp, expr[i], var[i])
#    end
#end
#
#function vipair(mp::MathPrgm, expr::Vector, var::JuMP.JuMPArray)
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        vipair(mp, expr[i], var[i])
#    end
#end

macro vipair(mp, expr, var)
    dump(expr, maxdepth=8)
    println(expr.head)
    if (expr.head == :vect || expr.head == :generator || expr.head == :comprehension)
        use_broadcast = true
    else
        use_broadcast = false
    end
    ast_sym2_esc!(expr)
    if use_broadcast
        dummyconstr = Expr(:call, :(.==), expr, 0)
    else
        dummyconstr = Expr(:call, :(==), expr, 0)
    end
    mmp = esc(mp)
    mvar = esc(var)
    return quote
        model = $(mmp).emp.model_ds
        cref = @constraint model $dummyconstr
        gidx = addequ!($mmp, cref)
        if (gidx isa Array)
            @assert length(gidx) == length($mvar)
            for i=1:length(gidx)
                $(mmp).matching[$(mvar)[i].col] = gidx[i]
            end
        else
            $(mmp).matching[$(esc(var)).col] = gidx
        end
        cref
    end
end

#macro vipair(mp, expr, var)
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        vipair(mp, expr[i], var[i])
#    end
#end
#
#macro vipair(mp, expr, var)
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        vipair(mp, expr[i], var[i])
#    end
#end
#
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

#function NLvipair(mp::MathPrgm, expr::Vector, var::Vector{JuMP.Variable})
#    @assert length(expr) == length(var)
#    for i in 1:length(var)
#        NLvipair(mp, expr[i], var[i])
#    end
#end


