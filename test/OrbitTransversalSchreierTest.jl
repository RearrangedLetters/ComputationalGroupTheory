using Test
using ComputationalGroupTheory

@testset "Orbit_1" begin
    σ = Permutation([2,3,4,1])
    @test orbit(σ, 1) == [1, 2, 3, 4]
    τ = Permutation([2,1])
    @test orbit(τ, 1) == [1, 2]
	S = [σ, τ]
	e = one(σ)
	Δ = orbit(S, e, *)
	@test length(Δ) == 24
    @test orbit(σ, 1) == [1, 2, 3, 4]
end