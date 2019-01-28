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

function ast_lookahead_loop_vars(ast_node)
    if (isa(ast_node, Array))
        return union(flatten(map(ast_lookahead_loop_vars, ast_node)))
    elseif (!isa(ast_node, Expr))
        return Vector{Symbol}(undef,0)
    end

    if (ast_node.head == :kw || ast_node.head == :generator || ast_node.head == :(=))
        return [ast_node.args[1]]
    end
    return Vector{Symbol}(undef,0)
end

function ast_sym2_esc!(ast_node, skipsym=Vector{Symbol}(undef,0))
    # TODO(xhub) move thus into a seperate function that would get called first
    if (isa(ast_node, Symbol))
        ast_node = esc(ast_node)
        return ast_node
    elseif (isa(ast_node, Array))
        return map(x -> ast_sym2_esc!(x, skipsym), ast_node)
    elseif (!isa(ast_node, Expr))
        return
    end

    # is we have those instances, do not quote
    if (ast_node.head == :call || ast_node.head == :callmacro || ast_node.head == :kw || ast_node.head == :(=))
        if (ast_node.head == :kw)
            push!(skipsym, ast_node.args[1])
        end
        s=2
    elseif (ast_node.head == :generator)
        # we have faith here that the iterator has a structure like
        # ast_node.args[2] ~ i=1:/N
        push!(skipsym, ast_node.args[2].args[1])
        # now visit ast_node.args[1] and ast_node.args[2].args[2] (i=1:N)
        if (typeof(ast_node.args[1]) != Expr)
            dump(ast_node)
            error("in ast_sym2_esc")
        end
        ast_sym2_esc!(ast_node.args[1], skipsym)
        if (typeof(ast_node.args[2].args[2]) != Expr)
            dump(ast_node.args[2])
            error("in ast_sym2_esc")
        end
        ast_sym2_esc!(ast_node.args[2].args[2], skipsym)
        return ast_node
    else
        s=1
    end

    for i=s:length(ast_node.args)
        if (typeof(ast_node.args[i]) == Symbol && !(ast_node.args[i] in skipsym))
            ast_node.args[i] = esc(ast_node.args[i])
        elseif (typeof(ast_node.args[i]) == Expr)
            ast_sym2_esc!(ast_node.args[i], skipsym)
        end
    end
    return ast_node
end

function ast_args_esc!(array_args)
    for i=1:length(array_args)
        if (typeof(array_args[i]) == Expr)
            array_args[i] = esc(array_args[i])
        end
    end
end

function ast_var_bracket_esc!(expr)
    if expr.head == :kw
        # We want to skip expr.args[1] since it is a loop variable
        ast_sym2_esc!(expr.args[2])
    else
        ast_sym2_esc!(expr)
    end
end

function esc_var!(ex)
    skipsym=Vector{Symbol}(undef,0)
    if (isa(ex, Symbol))
        varname = esc(ex)
    elseif (!isa(ex, Expr))
        error("Could not parse variable declaration $(ex); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
    elseif (ex.head == :ref || ex.head == :typed_vcat)
        # This is the case there we have var[...] or var[...; cond]
        # We can have var[i=1:3] or var[1:3]
        # - case var[1:3]: ex.args[2].head is call; ex.args[2].args is
        #   Array{Any}((3,))
        #     1: Symbol :
        #     2: Int64 1
        #     3: Int64 3
        #
        # case var[i=1:3]: ex.args[2].head is kw

        varname = esc(ex.args[1])
        # We want to skip escaping the inner loop var
        skipsym = ast_lookahead_loop_vars(ex.args[2:end])
        for i=2:length(ex.args)
            ex.args[i] = ast_sym2_esc!(ex.args[i], skipsym)
        end
    end
    return varname, skipsym
end

function esc_jump_id!(ex)
    if (typeof(ex) == Symbol)
        varname = esc(ex)
        ex = esc(ex)
        skipsym = Vector{Symbol}(undef,0)
    elseif (ex.head == :ref || ex.head == :typed_vcat)
        # This is the case there we have var[...]
        # We can have var[i=1:3] or var[1:3]
        # - case var[1:3]: ex.args[2].head is call; ex.args[2].args is
        #   Array{Any}((3,))
        #     1: Symbol :
        #     2: Int64 1
        #     3: Int64 3
        #
        # case var[i=1:3]: ex.args[2].head is kw

        varname, skipsym = esc_var!(ex)
        # We want to skip escaping the inner loop var
    elseif (ex.head == :call || ex.head == :comparison)
        # This is the case there we have var >= lb or ub >= var[...] >= lb
        if (ex.head == :comparison)
            posvar = 3
        else
            if (length(ex.args) == 3)
                posvar = 2
            else
                error("Could not parse variable declaration $(ex); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
            end
        end

        varname, skipsym = esc_var!(ex.args[posvar])
    else
        error("Could not parse variable declaration $(ex); please file a bug report at https://github.com/xhub/EMP.jl/issues/new")
    end

    return varname, skipsym
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
    # TODO(xhub) fix this with proper syntax
    # dummyconstr = Expr(:call, @variable, $(mp).emp.m, $(esc(args[1])))
    #then extend the arguments

    ar1 = args[1]
    varname, skipsym = esc_jump_id!(ar1)
    jump_call = :(JuMP.@variable mmp $args)

    if length(kwargs.args) > 0
        ast_sym2_esc!(kwargs, skipsym)
        append!(jump_call.args, kwargs.args)
    end

    # dump(args)
    return quote
        mmp = $(esc(mp)).emp.model_ds
        #varname = $(esc(args[1]))
        # TODO: use eval(Expr(:call, variable, ...))
        $varname = $jump_call
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
    ast_sym2_esc!(expr)
    dummyconstr = Expr(:call, :(==), expr, 0)
    mmp = esc(mp)
    return quote
        model = $(mmp).emp.model_ds
        cref = @constraint model $dummyconstr
        gidx = addequ!($mmp, cref)
        $(mmp).objequ = gidx
        $(mmp).sense = $(esc(sense))
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
    equname, skipsym = esc_jump_id!(name)
    ast_sym2_esc!(expr.args, skipsym)
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
    equname, skipsym = esc_jump_id!(name)
    ast_sym2_esc!(expr.args, skipsym)
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


