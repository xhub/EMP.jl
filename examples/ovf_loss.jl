# Data is from  ``Data analysis recipes: Fitting a model to data'' by Hogg, Bovy & Lang
# available on https://arxiv.org/abs/1008.4686.

if !isdef(:versions); versions = [1] end
if !isdef(:penalty_names); penalty_names = ["huber"] end
if !isdef(:ovf_formulations); ovf_formulations = ["equilibrium"] end

if VERSION >= v"0.7"
    using DelimitedFiles
end

N = length(y)

@testset "loss test: penalty = $penalty_name; version = $version; ovf_formulation = $ovf_formulation" for penalty_name in penalty_names, version in versions, ovf_formulation in ovf_formulations

solver = ReSHOPSolver("", Dict{String,Any}([("ovf_formulation", ovf_formulation)]))

m = JuMP.Model(solver=solver)

fit_model = EMP.Model(m)


if version == 1
    # primary unknonwn
    @variable(m, c1[j=1:1])
    @variable(m, d1)

    # OVF variable
    @variable(m, loss_var1)

    # OVF arguments
    @variable(m, fit1[i=1:N]);

    @constraint(m, fiteqn[i=1:N], (y[i] - sum(c1[j]*x[i, j] for j=1:M) - d1) / sigma_y[i] == fit1[i]);

    addovf!(fit_model, loss_var1, fit1, penalty_name, all_params[penalty_name])

    @objective(m, :Min, loss_var1)

    c = c1
    d = d1
    fit = fit1
    loss_var = loss_var1
elseif version == 2
    # primary unknonwn
    @variable(m, c2[j=1:1])
    @variable(m, d2)

    # OVF variable
    @variable(m, loss_var2)

    # OVF arguments
    @variable(m, fit2[i=1:N]);

    @NLconstraint(m, fiteqn[i=1:N], log(exp((y[i] - sum(c2[j]*x[i, j] for j=1:M) - d2) / sigma_y[i])) == fit2[i]);

    addovf!(fit_model, loss_var2, fit2, penalty_name, all_params[penalty_name])

    @objective(m, :Min, loss_var2)

    c = c2
    d = d2
    fit = fit2
    loss_var = loss_var2
elseif version == 3
    # primary unknonwn
    @variable(m, c3[j=1:1])
    @variable(m, d3)

    # OVF variable
    @variable(m, loss_var3)

    # OVF arguments
    @variable(m, fit3[i=1:N]);

    @constraint(m, fiteqn[i=1:N], fit3[i] == (y[i] - sum(c3[j]*x[i, j] for j=1:M) - d3) / sigma_y[i]);

    addovf!(fit_model, loss_var3, fit3, penalty_name, all_params[penalty_name])

    @objective(m, :Min, loss_var3)

    c = c3
    d = d3
    fit = fit3
    loss_var = loss_var3
elseif version == 4
    # primary unknonwn
    @variable(m, c4[j=1:1])
    @variable(m, d4)

    # OVF variable
    @variable(m, loss_var4)

    # OVF arguments
    @variable(m, fit4[i=1:N]);

    @NLconstraint(m, fiteqn[i=1:N], fit4[i] == log(exp((y[i] - sum(c4[j]*x[i, j] for j=1:M) - d4) / sigma_y[i])));

    addovf!(fit_model, loss_var4, fit4, penalty_name, all_params[penalty_name])

    @objective(m, :Min, loss_var4)

    c = c4
    d = d4
    fit = fit4
    loss_var = loss_var4
elseif version == 5
    # primary unknonwn
    @variable(m, c5[j=1:1])
    @variable(m, d5)

    # OVF variable
    @variable(m, loss_var5)

    # OVF arguments
    @variable(m, fit5[i=1:N]);

    @constraint(m, fiteqn[i=1:N], pi*(y[i] - sum(c5[j]*x[i, j] for j=1:M) - d5) / sigma_y[i] == pi*fit5[i]);

    addovf!(fit_model, loss_var5, fit5, penalty_name, all_params[penalty_name])

    @objective(m, :Min, loss_var5)

    c = c5
    d = d5
    fit = fit5
    loss_var = loss_var5
end

cmp_dir = joinpath(dirname(@__FILE__), "test_res")
unknown_ref = readdlm(joinpath(cmp_dir, "mcp_" * penalty_name * "_v1.out"))
file_ref = readdlm(joinpath(cmp_dir, "mcp_" * penalty_name * "_v1_fit.out"))

@test solve(m) == :Optimal

@test isapprox(getvalue(c), [unknown_ref[1]], rtol=1e-4)
@test isapprox(getvalue(d), unknown_ref[2], rtol=1e-4)
@test isapprox(getvalue(fit), file_ref[1:N], rtol=1e-4)

occursin("soft_hinge", penalty_name) && continue

@test isapprox(getvalue(loss_var), unknown_ref[3], rtol=1e-4)
#@test isapprox(getobjectivevalue(m), unknown_ref[3], rtol=1e-4)

end
