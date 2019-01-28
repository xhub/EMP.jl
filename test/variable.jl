#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for Variable
#############################################################################
using EMP, JuMP, Compat, Compat.Test, Compat.SparseArrays
import JuMP.repl, MathProgBase

if VERSION ≤ v"0.7-"
    dropdims(v; dims=nothing) = squeeze(v, dims)
end

@testset "Variables" begin
    @testset "constructors" begin
        # Constructors
        mcon = EMP.Model(); mp_mcon = MathPrgm(mcon)
        @variableMP(mp_mcon, nobounds)
        @variableMP(mp_mcon, lbonly >= 0)
        @variableMP(mp_mcon, ubonly <= 1)
        @variableMP(mp_mcon, 0 <= bothb <= 1)
        @variableMP(mp_mcon, 0 <= onerange[-5:5] <= 10)
        @variableMP(mp_mcon, onerangeub[-7:1] <= 10, Int)
        @variableMP(mp_mcon, manyrangelb[0:1,10:20,1:1] >= 2)
        @test getlowerbound(manyrangelb[0,15,1]) == 2
        s = ["Green","Blue"]
        @variableMP(mp_mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
        @test getupperbound(x[-4,"Green"]) == 5.5
        @test getvalue(x[-3,"Blue"]) == -2
        @test isequal(mcon.model_ds[:lbonly],lbonly)
        @test isequal(mcon.model_ds[:ubonly],ubonly)
        @test isequal(mcon.model_ds[:onerangeub][-7],onerangeub[-7])
        @variableMP(mp_mcon, lbonly)
        @test_throws ErrorException mcon.model_ds[:lbonly]
        @test_throws KeyError mcon.model_ds[:foo]
        # d exists and should be quoted ... idk how to deal with that
        # d = Dict()
        # @variableMP(mp_mcon, d["bar"][1:10] == 1)
        # @test getvalue(d["bar"][1]) == 1
        @test typeof(zero(nobounds)) == AffExpr
        @test typeof(one(nobounds)) == AffExpr
    end

    @testset "get and set bounds" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, 0 <= x <= 2)
        @test getlowerbound(x) == 0
        @test getupperbound(x) == 2
        setlowerbound(x, 1)
        @test getlowerbound(x) == 1
        setupperbound(x, 3)
        @test getupperbound(x) == 3
        @variableMP(mp_m, q, Bin)
        @test getlowerbound(q) == 0
        @test getupperbound(q) == 1
        @variableMP(mp_m, 0 <= y <= 1, Bin)
        @test getlowerbound(y) == 0
        @test getupperbound(y) == 1
        @variableMP(mp_m, fixedvar == 2)
        @test getvalue(fixedvar) == 2
        @test getlowerbound(fixedvar) == 2
        @test getupperbound(fixedvar) == 2
        JuMP.fix(fixedvar, 5)
        @test getvalue(fixedvar) == 5
        @test getlowerbound(fixedvar) == 5
        @test getupperbound(fixedvar) == 5
    end

    @testset "get and set values" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[1:3])
        x0 = collect(1:3)
        setvalue(x, x0)
        @test getvalue(x) == x0
        @test getvalue([x[1],x[2],x[3]]) == x0

        @variableMP(mp_m, y[1:3,1:2])
        @test_throws DimensionMismatch setvalue(y, collect(1:6))
    end

    @testset "get and set category" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[1:3])
        setcategory(x[2], :Int)
        @test getcategory(x[3]) == :Cont
        @test getcategory(x[2]) == :Int
    end

    # XXX this is smelly: can't really print it with JuMP ...
#    @testset "repeated elements in index set (issue #199)" begin
#        repeatmod = EMP.Model(); mp_repeatmod = MathPrgm(repeatmod)
#        s = [:x,:x,:y]
#        @variableMP(mp_repeatmod, x[s])
#        @test MathProgBase.numvar(repeatmod.model_ds) == 3
#    end

    @testset "condition in indexing" begin
        fa = repl[:for_all]
        inset, dots = repl[:in], repl[:dots]
        condmod = EMP.Model(); mp_condmod = MathPrgm(condmod)
        @variableMP(mp_condmod, x[i=1:10; iseven(i)])
        @variableMP(mp_condmod, y[j=1:10,k=3:2:9; isodd(j+k) && k <= 8])
        @test length(x.tupledict) == 5
        @test length(y.tupledict) == 15
        @test string(condmod.model_ds) == "Min 0\nSubject to\n x[i] $fa i $inset {1,2,$dots,9,10} s.t. iseven(i)\n y[j,k] $fa j $inset {1,2,$dots,9,10}, k $inset {3,5,7,9} s.t. isodd(j + k) and k <= 8\n"
    end

    @testset "@variable returning Array{Variable}" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[1:3,1:4,1:2])
        @variableMP(mp_m, y[1:0])
        @variableMP(mp_m, z[1:4])

        @test typeof(x) == Array{Variable,3}
        @test typeof(y) == Array{Variable,1}
        @test typeof(z) == Array{Variable,1}

        @test typeof(getvalue(x)) == Array{Float64,3}
        @test typeof(getvalue(y)) == Array{Float64,1}
        @test typeof(getvalue(z)) == Array{Float64,1}
    end

    @testset "getvalue on empty things" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[1:4,  1:0,1:3])   # Array{Variable}
        @variableMP(mp_m, y[1:4,  2:1,1:3]) # JuMPArray
        @variableMP(mp_m, z[1:4,Set(),1:3]) # JuMPDict

        @test getvalue(x) == Array{Float64}(undef, 4, 0, 3)
        @test typeof(getvalue(y)) <: JuMP.JuMPArray{Float64}
        @test JuMP.size(getvalue(y)) == (4,0,3)
        @test typeof(getvalue(z)) == JuMP.JuMPArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
        @test length(getvalue(z)) == 0
    end

# Slices three-dimensional JuMPContainer x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function sliceof(x, I, J, K)
    y = Array{Variable}(undef, length(I), length(J), length(K))

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

    @testset "Slices of JuMPArray (#684)" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[1:3, 1:4,1:2])
        @variableMP(mp_m, y[1:3,-1:2,3:4])
        @variableMP(mp_m, z[1:3,-1:2:4,3:4])
        @variableMP(mp_m, w[1:3,-1:2,[:red,"blue"]])

        @test x[:] == vec(sliceof(x, 1:3, 1:4, 1:2))
        @test x[:,:,:] == sliceof(x, 1:3, 1:4, 1:2)
        @test x[1,:,:] == sliceof(x, 1, 1:4, 1:2)
        @test x[1,:,2] == sliceof(x, 1, 1:4, 2)
        @test_throws BoundsError x[1,:,3]
        @test x[1:2,:,:] == sliceof(x, 1:2, 1:4, 1:2)
        @test x[1:2,:,2] == sliceof(x, 1:2, 1:4, 2)
        @test x[1:2,:,1:2] == sliceof(x, 1:2, 1:4, 1:2)
        @test_throws BoundsError x[1:2,:,1:3]

        @test y[:] == vec(sliceof(y, 1:3, -1:2, 3:4))
        @test y[:,:,:] == sliceof(y, 1:3, -1:2, 3:4)
        @test y[1,:,:] == sliceof(y, 1, -1:2, 3:4)
        @test y[1,:,4] == sliceof(y, 1, -1:2, 4)
        @test_throws ErrorException y[1,:,5]
        @test y[1:2,:,:] == sliceof(y, 1:2, -1:2, 3:4)
        @test y[1:2,:,4] == sliceof(y, 1:2, -1:2, 4)
        @test y[1:2,:,3:4] == sliceof(y, 1:2, -1:2, 3:4)
        @test_throws BoundsError y[1:2,:,1:3]

        @test z[:] == vec(sliceof(z, 1:3, -1:2:4, 3:4))
        @test z[:,1,:] == sliceof(z, 1:3, 1, 3:4)
        @test z[1,1,:] == sliceof(z, 1, 1, 3:4)
        @test_throws ErrorException z[:,5,3]
        @test z[1:2,1,:] == sliceof(z, 1:2, 1, 3:4)
        @test z[1:2,1,4] == sliceof(z, 1:2, 1, 4)
        @test z[1:2,1,3:4] == sliceof(z, 1:2, 1, 3:4)
        @test_throws BoundsError z[1:2,1,1:3]

        @test w[:] == vec(sliceof(w, 1:3, -1:2, [:red,"blue"]))
        @test_throws ErrorException w[:,:,:]
        @test w[1,:,"blue"] == sliceof(w, 1, -1:2, ["blue"])
        @test w[1,:,:red] == sliceof(w, 1, -1:2, [:red])
        @test_throws ErrorException w[1,:,"green"]
        @test w[1:2,:,"blue"] == sliceof(w, 1:2, -1:2, ["blue"])
        @test_throws ErrorException w[1:2,:,[:red,"blue"]]
    end

    @testset "Can't use end for indexing a JuMPContainer" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x[0:2,1:4])
        @variableMP(mp_m, y[i=1:4,j=1:4;true])
        @variableMP(mp_m, z[0:2])
        @test_throws ErrorException x[end,1]
        @test_throws ErrorException x[end-1]
        @test_throws ErrorException x[0,end-1]
        @test_throws MethodError y[end,end-1]
        @test_throws MethodError y[end,1]
        @test_throws ErrorException z[end]
    end

    @testset "Unsigned dimension lengths" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        t = UInt(4)
        @variableMP(mp_m, x[1:t])
        @constraintref(y[1:t])
        @test MathProgBase.numvar(m.model_ds) == 4
    end

    @testset "getvalue on sparse array (#889)" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        @variableMP(mp_m, x)
        @variableMP(mp_m, y)
        X = sparse([1, 3], [2, 3], [x, y])

        @test typeof(X) == SparseMatrixCSC{Variable, Int}

        setvalue(x, 1)
        setvalue(y, 2)

        Y = getvalue(X)
        @test typeof(Y) == SparseMatrixCSC{Float64, Int}
        @test Y == sparse([1, 3], [2, 3], [1, 2])
    end

    @testset "@constraintref does not work with JuMPArray #1329" begin
        m = EMP.Model(); mp_m = MathPrgm(m)
        cities = [:LA, :NY]
        @constraintref cref[cities]
        @test cref.indexsets[1] == cities
    end
end
