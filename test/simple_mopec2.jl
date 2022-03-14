using Test

@testset "simple mopec" begin

  include(joinpath(dirname(dirname(@__FILE__)), "examples", "simple_mopec.jl"))

  @test isapprox(value.(x)[:], [3,2,0], atol=1e-4)
  @test isapprox(value.(p)[:], [6,1,5], atol=1e-4)
  @test isapprox(value(y), 3, atol=1e-4)

#@test isapprox(objective_value(ag),5326.851310161077, atol=1e-5)

end
