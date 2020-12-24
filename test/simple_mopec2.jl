using JuMP, ReSHOP, EMP
using Test

@testset "simple mopec" begin

ATmat = [1., -1., -1.]
s = [.9, .1, 0]
b = [0, 5, 3]

mopec = EMPmaster()

ag = MathPrgm(mopec)
mkt = MathPrgm(mopec)

EquilibriumProblem(mopec, [ag, mkt])

n = 3

@variable(mkt, y >= 0)
@variable(ag, x[1:n] >= 0)
@variable(mkt, p[1:n] >= 0)

JuMP.fix(p[2], 1.; force=true)
JuMP.set_start_value.(x, 1.)

@NLobjective(ag, Max, sum(s[i] * log(x[i]) for i=1:n))
@constraint(ag, sum(p[i]*x[i] for i=1:n) <= sum(p[i]*b[i] for i=1:n) )

@constraint(mkt, [b[i] + ATmat[i]*y - x[i] for i=1:n] ⟂ p)
@constraint(mkt, sum(-ATmat[i]*p[i] for i=1:n) ⟂ y)

solveEMP(mopec)

@test isapprox(value.(x)[:], [3,2,0], atol=1e-4)
@test isapprox(value.(p)[:], [6,1,5], atol=1e-4)
@test isapprox(value(y), 3, atol=1e-4)

#@test isapprox(getobjectivevalue(m),5326.851310161077, atol=1e-5)

end
