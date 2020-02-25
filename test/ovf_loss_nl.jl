using JuMP, EMP, ReSHOP
using Test

# Data is from  ``Data analysis recipes: Fitting a model to data'' by Hogg, Bovy & Lang
# available on https://arxiv.org/abs/1008.4686.
y = [592 401 583 402 495 173 479 504 510 416 393 442 317 311 400 337 423 334 533 344]';
x = [201 244 47 287 203 58 210 202 198 158 165 201 157 131 166 160 186 125 218 146]';
sigma_y = [61 25 38 15 21 15 27 14 30 16 14 25 52 16 34 31 42 26 16 22]';

N = length(y)
M = 1

huber_params = Dict("kappa" => 1)
l1_params = Dict()
l2_params = Dict()
elastic_net_params = Dict("lambda" => 1)
hinge_params = Dict("epsilon" => 1)
vapnik_params = Dict("epsilon" => 1)
soft_hinge_params = Dict("epsilon" => 1, "kappa" => 1)
hubnik_params = Dict("epsilon" => 1, "kappa" => 1)

all_params = Dict(
"huber" => huber_params,
"l1" => l1_params,
"l2" => l2_params,
"elastic_net" => elastic_net_params,
"hinge" => hinge_params,
"vapnik" => vapnik_params,
"soft_hinge" => soft_hinge_params,
"hubnik" => hubnik_params
)

versions = 1:5
penalty_names = keys(all_params)
ovf_formulations = ["equilibrium", "dual"]


@testset "ovf_loss_nl" begin

#for penalty_name in penalty_names, version in versions

include(joinpath(dirname(dirname(@__FILE__)), "examples", "ovf_loss.jl"))

end
