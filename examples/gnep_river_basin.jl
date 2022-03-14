using JuMP, ReSHOP, EMP


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

ctx = EMPmaster()

ag = [MathPrgm(ctx) for i in 1:n]

EquilibriumProblem(ctx, ag)

# Add the variable to the MP
x = Vector{Any}(undef, n)
for i in 1:n
  var = @variable(ag[i], lower_bound=0)
  setindex!(x, var, i)
end

constr = Array{Any}(undef, n, 2)

for i in 1:n
    @objective(ag[i], Min, (c[1, i] + c[2, i]*x[i])*x[i] - (d1 - d2*sum(x[j] for j in 1:n))*x[i])
    for m in 1:2
        cons = @constraint(ag[i], sum(u[j, m]*ee[j]*x[j] for j in 1:3) <= K[m])
        setindex!(constr, cons, i, m)
    end
end

solveEMP(ctx)
