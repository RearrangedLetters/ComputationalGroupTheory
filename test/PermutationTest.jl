using Test

@testset "Permutation_1" begin
    σ = Permutation([2, 3, 4, 1])
    τ = Permutation([1, 3, 2])
    @test degree(σ) == 4
    @test inv(one(σ)) == one(σ)
	@test inv(σ) * σ == one(σ)
	@test τ * inv(τ) == one(τ)
	@test inv(σ * τ) == inv(τ) * inv(σ)
    σ = Permutation([2, 1, 3])
    @test σ * τ == Permutation([3, 1, 2])
end