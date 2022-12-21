using Test

@testset "Permutation_1" begin
    σ = Permutation([2, 3, 4, 1])
    τ = Permutation([1, 3, 2])
    @test degree(σ) == 4
    @test inv(one(σ)) == one(σ)
	@test inv(σ) * σ == one(σ)
	@test τ * inv(τ) == one(τ)
	@test inv(σ * τ) == inv(τ) * inv(σ)
    @test Permutation([2, 1, 3]) * τ == Permutation([3, 1, 2])
end

@testset "Permutation_2" begin
    σ = perm"(1, 2, 3, 4)(3, 4)"
    @test σ == Permutation([2, 3, 1, 4])
end

@testset "Permutation_3" begin
    σ = perm"(1, 2, 3)"
    @test σ^1 == σ
    @test σ^2 == σ * σ
    @test σ^2 == σ^σ
    @test σ^3 == one(σ)

    τ = perm"(1, 2)"
    η = perm"(2, 3)"
    @test τ * η == perm"(1, 3, 2)"
    @test η * τ == perm"(1, 2, 3)"
    @test τ^η == perm"(1, 3, 2)"
    @test η^τ == perm"(1, 2, 3)"
end