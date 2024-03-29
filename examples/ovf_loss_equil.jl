# Data is from  ``Data analysis recipes: Fitting a model to data'' by Hogg, Bovy & Lang
# available on https://arxiv.org/abs/1008.4686.

if !isdef(:penalty_names); penalty_names = ["huber"] end
if !isdef(:ovf_formulations); ovf_formulations = ["equilibrium"] end

@testset "Equilibrium loss test: penalty = $penalty_name; ovf_formulation = $ovf_formulation" for penalty_name in penalty_names, ovf_formulation in ovf_formulations

if ovf_formulation == "dual"
  if penalty_name == "hubnik"
    n_MP = 2
  else
    n_MP = 3
  end
else
  n_MP = 5
end

fit_model = EMPmaster()
m = get_JuMP_model(fit_model)
#fit_model.opts["ovf_formulation"] = ovf_formulation

mps = [MathPrgm(fit_model) for i=1:n_MP]

EquilibriumProblem(fit_model, mps)

# primary unknonwn
@variable(mps[1], c1[j=1:1])
@variable(mps[1], d1)

# OVF variable
@variable(mps[1], loss_var1)

# OVF arguments
@variable(mps[1], fit1[i=1:N]);

@constraint(mps[1], fiteqn1[i=1:N], (y[i] - sum(c1[j]*x[i, j] for j=1:M) - d1) / sigma_y[i] == fit1[i]);

addovf!(fit_model, loss_var1, fit1, penalty_name, all_params[penalty_name])

@objective(mps[1], Min, loss_var1)

# primary unknonwn
@variable(mps[2], c2[j=1:1])
@variable(mps[2], d2)

# OVF variable
@variable(mps[2], loss_var2)

# OVF arguments
@variable(mps[2], fit2[i=1:N]);

@NLconstraint(mps[2], fiteqn2[i=1:N], log(exp((y[i] - sum(c2[j]*x[i, j] for j=1:M) - d2) / sigma_y[i])) == fit2[i]);

addovf!(fit_model, loss_var2, fit2, penalty_name, all_params[penalty_name])

@objective(mps[2], Min, loss_var2)

if n_MP >= 3
# primary unknonwn
@variable(mps[3], c3[j=1:1])
@variable(mps[3], d3)

# OVF variable
@variable(mps[3], loss_var3)

# OVF arguments
@variable(mps[3], fit3[i=1:N]);

@constraint(mps[3], fiteqn3[i=1:N], fit3[i] == (y[i] - sum(c3[j]*x[i, j] for j=1:M) - d3) / sigma_y[i]);

addovf!(fit_model, loss_var3, fit3, penalty_name, all_params[penalty_name])

@objective(mps[3], Min, loss_var3)

elseif n_MP >= 4
  # primary unknonwn
  @variable(mps[4], c4[j=1:1])
  @variable(mps[4], d4)

  # OVF variable
  @variable(mps[4], loss_var4)

  # OVF arguments
  @variable(mps[4], fit4[i=1:N]);

  @NLconstraint(mps[4], fiteqn4[i=1:N], fit4[i] == log(exp((y[i] - sum(c4[j]*x[i, j] for j=1:M) - d4) / sigma_y[i])));

  addovf!(fit_model, loss_var4, fit4, penalty_name, all_params[penalty_name])

  @objective(mps[4], Min, loss_var4)

# primary unknonwn
  @variable(mps[5], c5[j=1:1])
  @variable(mps[5], d5)

  # OVF variable
  @variable(mps[5], loss_var5)

  # OVF arguments
  @variable(mps[5], fit5[i=1:N]);

  @constraint(mps[5], fiteqn5[i=1:N], pi*(y[i] - sum(c5[j]*x[i, j] for j=1:M) - d5) / sigma_y[i] == pi*fit5[i]);

  addovf!(fit_model, loss_var5, fit5, penalty_name, all_params[penalty_name])

  @objective(mps[5], Min, loss_var5)
end


cmp_dir = joinpath(dirname(@__FILE__), "test_res")
unknown_ref = readdlm(joinpath(cmp_dir, "mcp_" * penalty_name * "_v1.out"))
file_ref = readdlm(joinpath(cmp_dir, "mcp_" * penalty_name * "_v1_fit.out"))

solveEMP(fit_model)
@test EMP.termination_status(fit_model) == MOI.LOCALLY_SOLVED

occursin("hinge", penalty_name) && continue

@test isapprox(value.(c1), [unknown_ref[1]], rtol=1e-4)
@test isapprox(value(d1), unknown_ref[2], rtol=1e-4)
@test isapprox(value.(fit1), file_ref[1:N], rtol=1e-4)
@test isapprox(value(loss_var1), unknown_ref[3], rtol=1e-4)

@test isapprox(value.(c2), [unknown_ref[1]], rtol=1e-4)
@test isapprox(value(d2), unknown_ref[2], rtol=1e-4)
@test isapprox(value.(fit2), file_ref[1:N], rtol=1e-4)
@test isapprox(value(loss_var2), unknown_ref[3], rtol=1e-4)

if n_MP >= 3
@test isapprox(value.(c3), [unknown_ref[1]], rtol=1e-4)
@test isapprox(value(d3), unknown_ref[2], rtol=1e-4)
@test isapprox(value.(fit3), file_ref[1:N], rtol=1e-4)
@test isapprox(value(loss_var3), unknown_ref[3], rtol=1e-4)

elseif n_MP >= 4
  @test isapprox(value.(c4), [unknown_ref[1]], rtol=1e-4)
  @test isapprox(value(d4), unknown_ref[2], rtol=1e-4)
  @test isapprox(value.(fit4), file_ref[1:N], rtol=1e-4)
  @test isapprox(value(loss_var4), unknown_ref[3], rtol=1e-4)

  @test isapprox(value.(c5), [unknown_ref[1]], rtol=1e-4)
  @test isapprox(value(d5), unknown_ref[2], rtol=1e-4)
  @test isapprox(value.(fit5), file_ref[1:N], rtol=1e-4)
  @test isapprox(value(loss_var5), unknown_ref[3], rtol=1e-4)
end

end
