using JuMP, JAMSDWriter, EMP
using Base.Test

@testset "gnep river basin" begin

K = [100 100]
d1 = 3
d2 = 0.01
ee = [0.5, 0.25, 0.75]

c = [0.1 0.12 0.15;
    0.01 0.05 0.01]

u = [6.5    4.583;
     5.0    6.250;
     5.5    3.750]

n = 3

ctx = EMP.Model()

jump_model = ctx.model_ds

ag = [MathPrgm(ctx) for i in 1:n]

@variable(jump_model, x[i=1:n] >= 0)

# Add the variable to the MP
for i in 1:n
    addvar!(ag[i], x[i])
end

constr = Array{Any}(n, 2)

for i in 1:n
    @objectiveMP(ag[i], :Min, (c[1, i] + c[2, i]*x[i])*x[i] - (d1 - d2*sum(x[j] for j in 1:n))*x[i])
    for m in 1:2
        cons = @constraintMP(ag[i], sum(u[j, m]*ee[j]*x[j] for j in 1:3) <= K[m])
        setindex!(constr, cons, i, m)
    end
end

@test solve(jump_model) == :Optimal

@test isapprox(getvalue(x), [0., 6.47333, 22.2808], atol=1e-4)
@test isapprox(getdual(constr[1,:]), [-0.8038, 0], atol=1e-3)
@test isapprox(getdual(constr[2,:]), [-1.5043, 0], atol=1e-3)
@test isapprox(getdual(constr[3,:]), [-0.4592, 0], atol=1e-3)

end
