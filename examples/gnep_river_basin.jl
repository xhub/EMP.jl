using JuMP, JAMSDWriter, EMP

@testset "river basin (GNEP)" begin

K = [100 100]
d1 = 3
d2 = 0.01
ee = [0.5, 0.25, 0.75]

c = [0.1 0.12 0.15;
    0.01 0.05 0.01]

u = [6.5    4.583;
     5.0    6.250;
     5.5    3.750]

solver = JAMSDWriter.JAMSDSolver()

n = 3

jump_model = JuMP.Model(solver=solver)

ctx = EMP.Model(jump_model)

ag = [EMP.MathPrgm(ctx) for i in 1:n]

@variable(jump_model, x[i=1:n] >= 0)
for i in 1:n
    EMP.addvar!(ag[i], x[i])
end

for i in 1:n
    objectiveMP(ag[i], :Min, (c[1, i] + c[2, i]*x[i])*x[i] - (d1 - d2*sum(x[j] for j in 1:n))*x[i])
    for m in 1:2
        @constraintMP(ag[i], sum(u[j, m]*ee[j]*x[j] for j in 1:3) <= K[m])
    end
end

solveEMP(ctx)

println("$(jump_model.internalModel.inner.solution)")

end
