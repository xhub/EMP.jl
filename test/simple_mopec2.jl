using JuMP, ReSHOP, EMP
using Compat.Test

@testset "simple mopec" begin

ATmat = [1., -1., -1.]
s = [.9, .1, 0]
b = [0, 5, 3]

mopec = EMP.Model()

ag = MathPrgm(mopec)
mkt = MathPrgm(mopec)

EquilibriumProblem(mopec, [ag, mkt])

n = 3

@variableMP(mkt, y >= 0)
@variableMP(ag, x[1:n] >= 0)
@variableMP(mkt, p[1:n] >= 0)

JuMP.fix(p[2], 1.)
JuMP.setvalue(x[1:n], ones(n))

@NLobjectiveMP(ag, :Max, sum(s[i] * log(x[i]) for i=1:n))
@constraintMP(ag, sum(p[i]*x[i] for i=1:n) <= sum(p[i]*b[i] for i=1:n) )

vipair(mkt, [b[i] + ATmat[i]*y - x[i] for i=1:n], p)
vipair(mkt, sum(-ATmat[i]*p[i] for i=1:n), y)

@test solveEMP(mopec) == :Optimal

@test isapprox(getvalue(x), [3,2,0], atol=1e-4)
@test isapprox(getvalue(p), [6,1,5], atol=1e-4)
@test isapprox(getvalue(y), 3, atol=1e-4)

#@test isapprox(getobjectivevalue(m),5326.851310161077, atol=1e-5)

end
