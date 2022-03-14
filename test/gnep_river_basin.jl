using Test

@testset "gnep river basin" begin

  include(joinpath(dirname(dirname(@__FILE__)), "examples", "gnep_river_basin.jl"))

  @test isapprox(value.(x), [0., 6.47333, 22.2808], atol=1e-4)
  @test isapprox(dual.(constr[1,:]), [-0.8038, 0], atol=1e-3)
  @test isapprox(dual.(constr[2,:]), [-1.5043, 0], atol=1e-3)
  @test isapprox(dual.(constr[3,:]), [-0.4592, 0], atol=1e-3)

end
