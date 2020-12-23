#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/jump-dev/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for VariableRef
#############################################################################

using JuMP
using EMP

import LinearAlgebra: Symmetric
using Test

include(joinpath(@__DIR__, "utilities.jl"))
#@static if !(:JuMPExtension in names(Main))
#    include("JuMPExtension.jl")
#end

VariableRefType = JuMP.VariableRef

function test_variable_name(variable, name)
    @test name == @inferred JuMP.name(variable)
    @test variable == JuMP.variable_by_name(JuMP.owner_model(variable), name)
end

# Slices three-dimensional DenseAxisArray x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function _sliceof_util(VariableRefType, x, I, J, K)
    y = Array{VariableRefType}(undef, length(I), length(J), length(K))

    ii = 1
    jj = 1
    kk = 1
    for i in I
        for j in J
            for k in K
                y[ii,jj,kk] = x[i,j,k]
                kk += 1
            end
            jj += 1
            kk = 1
        end
        ii += 1
        jj = 1
    end
    idx = [length(I)==1, length(J)==1, length(K)==1]
    dropdims(y, dims=tuple(findall(idx)...))
end

function test_variable_no_bound(ModelType, VariableRefType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, nobounds)
    @test !JuMP.has_lower_bound(nobounds)
    @test !JuMP.has_upper_bound(nobounds)
    @test !JuMP.is_fixed(nobounds)
    test_variable_name(nobounds, "nobounds")
    @test zero(nobounds) isa JuMP.GenericAffExpr{Float64, VariableRefType}
    @test one(nobounds) isa JuMP.GenericAffExpr{Float64, VariableRefType}
end

function test_variable_lower_bound_rhs(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, lbonly >= 0, Bin)
    @test JuMP.has_lower_bound(lbonly)
    @test 0.0 == @inferred JuMP.lower_bound(lbonly)
    @test !JuMP.is_fixed(lbonly)
    @test JuMP.is_binary(lbonly)
    @test !JuMP.is_integer(lbonly)
    @test isequal(mp[:lbonly], lbonly)
    JuMP.delete_lower_bound(lbonly)
    @test !JuMP.has_lower_bound(lbonly)
    # Name already used
    @test_throws ErrorException @variable(mp, lbonly)
end

function test_variable_lower_bound_lhs(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, 0 <= lblhs, Bin)
    @test JuMP.has_lower_bound(lblhs)
    @test 0.0 == @inferred JuMP.lower_bound(lblhs)
    @test !JuMP.is_fixed(lblhs)
    @test JuMP.is_binary(lblhs)
    @test !JuMP.is_integer(lblhs)
    @test isequal(mp[:lblhs], lblhs)
end

function test_variable_upper_bound_rhs(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, ubonly <= 1, Int)
    @test !JuMP.has_lower_bound(ubonly)
    @test JuMP.has_upper_bound(ubonly)
    @test 1.0 == @inferred JuMP.upper_bound(ubonly)
    @test !JuMP.is_fixed(ubonly)
    @test !JuMP.is_binary(ubonly)
    @test JuMP.is_integer(ubonly)
    @test isequal(mp[:ubonly], ubonly)
    JuMP.delete_upper_bound(ubonly)
    @test !JuMP.has_upper_bound(ubonly)
end

function test_variable_upper_bound_lhs(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, 1 >= ublhs, Int)
    @test !JuMP.has_lower_bound(ublhs)
    @test JuMP.has_upper_bound(ublhs)
    @test 1.0 == @inferred JuMP.upper_bound(ublhs)
    @test !JuMP.is_fixed(ublhs)
    @test !JuMP.is_binary(ublhs)
    @test JuMP.is_integer(ublhs)
    @test isequal(mp[:ublhs],ublhs)
end

function test_variable_interval(ModelType)
    function has_bounds(var, lb, ub)
        @test JuMP.has_lower_bound(var)
        @test lb == @inferred JuMP.lower_bound(var)
        @test JuMP.has_upper_bound(var)
        @test ub == @inferred JuMP.upper_bound(var)
        @test !JuMP.is_fixed(var)
    end
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, 0 <= bothb1 <= 1)
    has_bounds(bothb1, 0.0, 1.0)
    @variable(mp, 0 ≤  bothb2 ≤  1)
    has_bounds(bothb2, 0.0, 1.0)
    @variable(mp, 1 >= bothb3 >= 0)
    has_bounds(bothb3, 0.0, 1.0)
    @variable(mp, 1 ≥  bothb4 ≥  0)
    has_bounds(bothb4, 0.0, 1.0)
    @test_macro_throws ErrorException @variable(mp, 1 ≥ bothb5 ≤ 0)
    @test_macro_throws ErrorException @variable(mp, 1 ≤ bothb6 ≥ 0)
end

function test_variable_fix(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, fixed == 1.0)
    @test !JuMP.has_lower_bound(fixed)
    @test !JuMP.has_upper_bound(fixed)
    @test JuMP.is_fixed(fixed)
    @test 1.0 == @inferred JuMP.fix_value(fixed)
    JuMP.unfix(fixed)
    @test !JuMP.is_fixed(fixed)
    JuMP.set_lower_bound(fixed, 0.0)
    @test_throws Exception JuMP.fix(fixed, 1.0)
    JuMP.fix(fixed, 1.0; force = true)
    @test !JuMP.has_lower_bound(fixed)
    @test !JuMP.has_upper_bound(fixed)
    @test JuMP.is_fixed(fixed)
    @test 1.0 == @inferred JuMP.fix_value(fixed)
end

function test_variable_custom_index_sets(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, onerangeub[-7:1] <= 10, Int)
    @variable(mp, manyrangelb[0:1, 10:20, 1:1] >= 2)
    @test JuMP.has_lower_bound(manyrangelb[0, 15, 1])
    @test 2 == @inferred JuMP.lower_bound(manyrangelb[0, 15, 1])
    @test !JuMP.has_upper_bound(manyrangelb[0, 15, 1])

    s = ["Green","Blue"]
    @variable(mp, x[i=-10:10, s] <= 5.5, Int, start=i+1)
    @test 5 == @inferred JuMP.upper_bound(x[-4, "Green"])
    test_variable_name(x[-10, "Green"], "x[-10,Green]")
    @test JuMP.start_value(x[-3, "Blue"]) == -2
    @test isequal(mp[:onerangeub][-7], onerangeub[-7])
    @test_throws KeyError mp[:foo]
end

function test_variable_anonymous(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @test_throws ErrorException @variable(mp, [(0, 0)])  # #922
    x = @variable(mp, [(0, 2)])
    @test "" == @inferred JuMP.name(x[0])
    @test "" == @inferred JuMP.name(x[2])
end

function test_variable_is_valid_delete(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x)
    @test JuMP.is_valid(mp, x)
    JuMP.delete(mp, x)
    @test !JuMP.is_valid(mp, x)
    second_m = EMPmaster(); mp = MathPrgm(m)
    @test_throws Exception JuMP.delete(second_mp, x)
end

function test_variable_bounds_set_get(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, 0 <= x <= 2)
    @test 0 == @inferred JuMP.lower_bound(x)
    @test 2 == @inferred JuMP.upper_bound(x)
    set_lower_bound(x, 1)
    @test 1 == @inferred JuMP.lower_bound(x)
    set_upper_bound(x, 3)
    @test 3 == @inferred JuMP.upper_bound(x)
    @variable(mp, q, Bin)

    @variable(mp, 0 <= y <= 1, Bin)
    @test 0 == @inferred JuMP.lower_bound(y)
    @test 1 == @inferred JuMP.upper_bound(y)

    @variable(mp, fixedvar == 2)
    @test 2.0 == @inferred JuMP.fix_value(fixedvar)
    JuMP.fix(fixedvar, 5)
    @test 5 == @inferred JuMP.fix_value(fixedvar)
    @test_throws Exception JuMP.lower_bound(fixedvar)
    @test_throws Exception JuMP.upper_bound(fixedvar)
end

function test_variable_starts_set_get(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:3])
    x0 = collect(1:3)
    JuMP.set_start_value.(x, x0)
    @test JuMP.start_value.(x) == x0
    @test JuMP.start_value.([x[1],x[2],x[3]]) == x0

    @variable(mp, y[1:3,1:2])
    @test_throws DimensionMismatch JuMP.set_start_value.(y, collect(1:6))
end

function test_variable_integrality_set_get(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:3])

    JuMP.set_integer(x[2])
    JuMP.set_integer(x[2])  # test duplicated call
    @test JuMP.is_integer(x[2])
    JuMP.unset_integer(x[2])
    @test !JuMP.is_integer(x[2])

    JuMP.set_binary(x[1])
    JuMP.set_binary(x[1])  # test duplicated call
    @test JuMP.is_binary(x[1])
    @test_throws Exception JuMP.set_integer(x[1])
    JuMP.unset_binary(x[1])
    @test !JuMP.is_binary(x[1])

    @variable(mp, y, binary = true)
    @test JuMP.is_binary(y)
    @test_throws Exception JuMP.set_integer(y)
    JuMP.unset_binary(y)
    @test !JuMP.is_binary(y)

    @variable(mp, z, integer = true)
    @test JuMP.is_integer(z)
    @test_throws Exception JuMP.set_binary(z)
    JuMP.unset_integer(z)
    @test !JuMP.is_integer(z)
end

function test_variable_repeated_elements(ModelType)
    # Tests repeated elements in index set throw error (JuMP issue #199).
    m = EMPmaster(); mp = MathPrgm(m)
    index_set = [:x,:x,:y]
    @test_throws ErrorException (
        @variable(mp, unused_variable[index_set], container=DenseAxisArray))
    @test_throws ErrorException (
        @variable(mp, unused_variable[index_set], container=SparseAxisArray))
    @test_throws ErrorException (
        @variable(mp, unused_variable[index_set, [1]], container=DenseAxisArray))
    @test_throws ErrorException (
        @variable(mp, unused_variable[index_set, [1]], container=SparseAxisArray))
end

function test_variable_oneto_index_set(ModelType, VariableRefType)
    # Tests that Base.OneTo can be used in index set (JuMP issue #933).
    m = EMPmaster(); mp = MathPrgm(m)
    auto_var = @variable(mp, [Base.OneTo(3), 1:2], container=Auto)
    @test auto_var isa Matrix{VariableRefType}
    @test (3, 2) == @inferred size(auto_var)
    array_var = @variable(mp, [Base.OneTo(3), 1:2], container=Array)
    @test array_var isa Matrix{VariableRefType}
    @test (3, 2) == @inferred size(array_var)
    denseaxisarray_var = @variable(mp, [Base.OneTo(3), 1:2], container=DenseAxisArray)
    @test denseaxisarray_var isa JuMP.Containers.DenseAxisArray{VariableRefType}
    @test length.(axes(denseaxisarray_var)) == (3, 2)
end

function test_variable_base_name_in_macro(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, normal_var)
    test_variable_name(normal_var, "normal_var")
    no_indices = @variable(mp, base_name="foo")
    test_variable_name(no_indices, "foo")
    # Note that `z` will be ignored in name.
    indices = @variable(mp, z[i=2:3], base_name="t")
    test_variable_name(indices[2], "t[2]")
    test_variable_name(indices[3], "t[3]")
end

function test_variable_name(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x)
    test_variable_name(x, "x")
    JuMP.set_name(x, "y")
    @test JuMP.variable_by_name(mp, "x") isa Nothing
    test_variable_name(x, "y")
    y = @variable(mp, base_name="y")
    err(name) = ErrorException("Multiple variables have the name $name.")
    @test_throws err("y") JuMP.variable_by_name(mp, "y")
    JuMP.set_name(y, "x")
    test_variable_name(x, "y")
    test_variable_name(y, "x")
    JuMP.set_name(x, "x")
    @test_throws err("x") JuMP.variable_by_name(mp, "x")
    @test JuMP.variable_by_name(mp, "y") isa Nothing
    JuMP.set_name(y, "y")
    test_variable_name(x, "x")
    test_variable_name(y, "y")
end

function test_variable_condition_in_indexing(ModelType)
    function test_one_dim(x)
        @test 5 == @inferred length(x)
        for i in 1:10
            if iseven(i)
                @test haskey(x, i)
            else
                @test !haskey(x, i)
            end
        end
    end

    function test_two_dim(y)
        @test 15 == @inferred length(y)
        for j in 1:10, k in 3:2:9
            if isodd(j+k) && k <= 8
                @test haskey(y, (j,k))
            else
                @test !haskey(y, (j,k))
            end
        end
    end

    m = EMPmaster(); mp = MathPrgm(m)
    # Parses as ref on 0.7.
    @variable(mp, named_one_dim[i=1:10; iseven(i)])
    test_one_dim(named_one_dim)
    # Parses as vcat on 0.7.
    anon_one_dim = @variable(mp, [i=1:10; iseven(i)])
    test_one_dim(anon_one_dim)
    # Parses as typed_vcat on 0.7.
    @variable(mp, named_two_dim[j=1:10, k=3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(named_two_dim)
    # Parses as vect on 0.7.
    anon_two_dim = @variable(mp, [j=1:10, k=3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(anon_two_dim)
end

function test_variable_macro_return_type(ModelType, VariableRefType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:3, 1:4, 1:2], start=0.0)
    @test typeof(x) == Array{VariableRefType,3}
    @test typeof(JuMP.start_value.(x)) == Array{Float64,3}

    @variable(mp, y[1:0], start=0.0)
    @test typeof(y) == Vector{VariableRefType}
    # No type to infer for an empty collection.
    @test typeof(JuMP.start_value.(y)) == Vector{Union{Nothing, Float64}}

    @variable(mp, z[1:4], start = 0.0)
    @test typeof(z) == Vector{VariableRefType}
    @test typeof(JuMP.start_value.(z)) == Vector{Float64}
end

function test_variable_start_value_on_empty(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:4,  1:0,1:3], start = 0)  # Array{VariableRef}
    @variable(mp, y[1:4,  2:1,1:3], start = 0)  # DenseAxisArray
    @variable(mp, z[1:4,Set(),1:3], start = 0)  # SparseAxisArray

    @test JuMP.start_value.(x) == Array{Float64}(undef, 4, 0, 3)
    # TODO: Decide what to do here. I don't know if we still need to test this
    #       given broadcast syntax.
    # @test typeof(JuMP.start_value(y)) <: JuMP.DenseAxisArray{Float64}
    # @test JuMP.size(JuMP.start_value(y)) == (4,0,3)
    # @test typeof(JuMP.start_value(z)) ==
    #   JuMP.DenseAxisArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
    # @test length(JuMP.start_value(z)) == 0
end

function test_variable_denseaxisarray_slices(ModelType, VariableRefType)
    # Test slicing DenseAxisArrays (JuMP issue #684).
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:3, 1:4, 1:2], container=DenseAxisArray)
    @variable(mp, y[1:3, -1:2, 3:4])
    @variable(mp, z[1:3, -1:2:4, 3:4])
    @variable(mp, w[1:3, -1:2,[:red, "blue"]])

    #@test x[:] == vec(_sliceof_util(VariableRefType, x, 1:3, 1:4, 1:2))
    @test x isa JuMP.Containers.DenseAxisArray
    @test x[:, :, :].data == _sliceof_util(VariableRefType, x, 1:3, 1:4, 1:2)
    @test x[1, :, :].data == _sliceof_util(VariableRefType, x, 1, 1:4, 1:2)
    @test x[1, :, 2].data == _sliceof_util(VariableRefType, x, 1, 1:4, 2)
    @test_throws KeyError x[1, :, 3]
    #@test x[1:2,:,:].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 1:2)
    #@test x[1:2,:,2].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 2)
    #@test x[1:2,:,1:2].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 1:2)
    @test_throws KeyError x[1:2, :, 1:3]

    #@test y[:] == vec(_sliceof_util(VariableRefType, y, 1:3, -1:2, 3:4))
    @test y[:, :, :].data == _sliceof_util(VariableRefType, y, 1:3, -1:2, 3:4)
    @test y[1, :, :].data == _sliceof_util(VariableRefType, y, 1, -1:2, 3:4)
    @test y[1, :, 4].data == _sliceof_util(VariableRefType, y, 1, -1:2, 4)
    @test_throws KeyError y[1, :, 5]
    # @test y[1:2,:,:] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test y[1:2,:,4] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 4)
    # @test y[1:2,:,3:4] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test_throws BoundsError y[1:2,:,1:3]

    #@test z[:] == vec(_sliceof_util(VariableRefType, z, 1:3, -1:2:4, 3:4))
    @test z[:, 1, :].data == _sliceof_util(VariableRefType, z, 1:3, 1, 3:4)
    @test z[1, 1, :].data == _sliceof_util(VariableRefType, z, 1, 1, 3:4)
    @test_throws KeyError z[:, 5, 3]
    # @test z[1:2,1,:] == _sliceof_util(VariableRefType, z, 1:2, 1, 3:4)
    # @test z[1:2,1,4] == _sliceof_util(VariableRefType, z, 1:2, 1, 4)
    # @test z[1:2,1,3:4] == _sliceof_util(VariableRefType, z, 1:2, 1, 3:4)
    # @test_throws BoundsError z[1:2,1,1:3]

    #@test w[:] == vec(_sliceof_util(VariableRefType, w, 1:3, -1:2, [:red,"blue"]))
    @test w[:, :, :] == w
    @test w[1, :, "blue"].data == _sliceof_util(VariableRefType, w, 1, -1:2, ["blue"])
    @test w[1, :, :red].data == _sliceof_util(VariableRefType, w, 1, -1:2, [:red])
    @test_throws KeyError w[1, :, "green"]
    # @test w[1:2,:,"blue"] == _sliceof_util(VariableRefType, w, 1:2, -1:2, ["blue"])
    # @test_throws ErrorException w[1:2,:,[:red,"blue"]]
end

function test_variable_end_indexing(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[0:2, 1:4])
    @variable(mp, z[0:2])
    @test x[end,1] == x[2, 1]
    @test x[0, end-1] == x[0, 3]
    @test z[end] == z[2]
    # TODO: It is redirected to x[11] as it is the 11th element but linear
    #       indexing is not supported
    @test_throws BoundsError x[end-1]
end

function test_variable_unsigned_index(ModelType)
    # Tests unsigned int can be used to construct index set (JuMP issue #857).
    m = EMPmaster(); mp = MathPrgm(m)
    t = UInt(4)
    @variable(mp, x[1:t])
    @test 4 == @inferred num_variables(mp)
end

function test_variable_symmetric(ModelType)
    m = EMPmaster(); mp = MathPrgm(m)

    @variable mp x[1:2, 1:2] Symmetric
    @test x isa Symmetric
    @test x[1, 2] === x[2, 1]
    @test mp[:x] === x

    y = @variable mp [1:2, 1:2] Symmetric
    @test y isa Symmetric
    @test y[1, 2] === y[2, 1]
end

function variables_test(ModelType::Type{<:JuMP.AbstractModel},
                        VariableRefType::Type{<:JuMP.AbstractVariableRef})
    @testset "Variable name" begin
        test_variable_name(ModelType)
    end

    @testset "Constructors" begin
        test_variable_no_bound(ModelType, VariableRefType)
        test_variable_lower_bound_rhs(ModelType)
        test_variable_lower_bound_lhs(ModelType)
        test_variable_upper_bound_rhs(ModelType)
        test_variable_upper_bound_lhs(ModelType)
        test_variable_interval(ModelType)
        test_variable_fix(ModelType)
        test_variable_custom_index_sets(ModelType)
        test_variable_anonymous(ModelType)
    end

    @testset "isvalid and delete variable" begin
        test_variable_is_valid_delete(ModelType)
    end

    @testset "get and set bounds" begin
        test_variable_bounds_set_get(ModelType)
    end

    @testset "get and set start" begin
        test_variable_starts_set_get(ModelType)
    end

    @testset "get and set integer/binary" begin
        test_variable_integrality_set_get(ModelType)
    end

    @testset "repeated elements in index set (issue #199)" begin
        test_variable_repeated_elements(ModelType)
    end

    @testset "Base.OneTo as index set (#933)" begin
        test_variable_oneto_index_set(ModelType, VariableRefType)
    end

    @testset "base_name= in @variable" begin
        test_variable_base_name_in_macro(ModelType)
    end

    @testset "condition in indexing" begin
        test_variable_condition_in_indexing(ModelType)
    end

    @testset "@variable returning Array{VariableRef}" begin
        test_variable_macro_return_type(ModelType, VariableRefType)
    end

    @testset "start_value on empty things" begin
        test_variable_start_value_on_empty(ModelType)
    end

    @testset "Slices of DenseAxisArray (#684)" begin
        test_variable_denseaxisarray_slices(ModelType, VariableRefType)
    end

    @testset "end for indexing a DenseAxisArray" begin
        test_variable_end_indexing(ModelType)
    end

    @testset "Unsigned dimension lengths (#857)" begin
        test_variable_unsigned_index(ModelType)
    end

end

@testset "Variables for JuMP.Model" begin
    variables_test(Model, VariableRef)
    @testset "all_variables" begin
        m = EMPmaster(); mp = MathPrgm(m)
        @variable(mp, x)
        @variable(mp, y)
        @test [x, y] == @inferred JuMP.all_variables(mp)
    end
#    @testset "@variables" begin
#        m = EMPmaster(); mp = MathPrgm(m)
#        @variables model begin
#            0 ≤ x[i=1:2] ≤ i
#            y ≥ 2, Int, (start = 0.7)
#            z ≤ 3, (start = 10)
#            q, (Bin, start = 0.5)
#        end
#
#        @test "x[1]" == @inferred JuMP.name(x[1])
#        @test 0 == @inferred JuMP.lower_bound(x[1])
#        @test 1 == @inferred JuMP.upper_bound(x[1])
#        @test !JuMP.is_binary(x[1])
#        @test !JuMP.is_integer(x[1])
#        @test JuMP.start_value(x[1]) === nothing
#
#        @test "x[2]" == @inferred JuMP.name(x[2])
#        @test 0 == @inferred JuMP.lower_bound(x[2])
#        @test 2 == @inferred JuMP.upper_bound(x[2])
#        @test !JuMP.is_binary(x[2])
#        @test !JuMP.is_integer(x[2])
#        @test JuMP.start_value(x[2]) === nothing
#
#        @test "y" == @inferred JuMP.name(y)
#        @test 2 == @inferred JuMP.lower_bound(y)
#        @test !JuMP.has_upper_bound(y)
#        @test !JuMP.is_binary(y)
#        @test JuMP.is_integer(y)
#        @test JuMP.start_value(y) === 0.7
#
#        @test "z" == @inferred JuMP.name(z)
#        @test !JuMP.has_lower_bound(z)
#        @test 3 == @inferred JuMP.upper_bound(z)
#        @test !JuMP.is_binary(z)
#        @test !JuMP.is_integer(z)
#        @test JuMP.start_value(z) === 10.0
#
#        @test "q" == @inferred JuMP.name(q)
#        @test !JuMP.has_lower_bound(q)
#        @test !JuMP.has_upper_bound(q)
#        @test JuMP.is_binary(q)
#        @test !JuMP.is_integer(q)
#        @test JuMP.start_value(q) === 0.5
#    end
end

#@testset "Variables for JuMPExtension.MyModel" begin
#    variables_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
#end

@testset "Dual from Variable" begin
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x == 0)
    exception = ErrorException(
        "To query the dual variables associated with a variable bound, first " *
        "obtain a constraint reference using one of `UpperBoundRef`, `LowerBoundRef`, " *
        "or `FixRef`, and then call `dual` on the returned constraint reference.\nFor " *
        "example, if `x <= 1`, instead of `dual(x)`, call `dual(UpperBoundRef(x))`.")
    @test_throws exception JuMP.dual(x)
end

@testset "value on containers" begin
    m = EMPmaster(); mp = MathPrgm(m)
    @variable(mp, x[1:2])
    exception = ErrorException(
        "`JuMP.value` is not defined for collections of JuMP types. Use " *
        "Julia's broadcast syntax instead: `JuMP.value.(x)`.")
    @test_throws exception JuMP.value(x)
end
