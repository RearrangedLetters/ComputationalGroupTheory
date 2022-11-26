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

@testset "Orbit_2" begin
    σ = perm"(1, 2, 3, 4)"
    @test orbit(σ, 1) == [1, 2, 3, 4]
    τ = perm"(1, 2)"
    @test orbit(τ, 1) == [1, 2]
	S = [σ, τ]
	e = one(σ)
	Δ = orbit(S, e, *)
	@test length(Δ) == 24
    @test orbit(σ, 1) == [1, 2, 3, 4]
end

@testset "FactoredTransversal_1" begin
    σ = Permutation([1, 3, 4, 2])
    τ = Permutation([1, 2, 4, 5, 3])
    x = 2
	Δ, T = transversalFactored([σ, τ], x)
	@test length(Δ) == 4
	for δ in Δ
		@test x^prod(T[δ]) == δ
	end
	Δ, T
end

@testset "FactoredTransversal_2" begin
    σ = Permutation([1,4,2,3])
    τ = Permutation([2,3,1])
    x = one(σ)
    S = [σ, τ]
    e = one(σ)
    Δ, T = transversalFactored(S, e, *)
    @test length(Δ) == 12
    for g in Δ
        @test g == prod(T[g])
    end
end

@testset "SchreierTransversal_1" begin
    σ = Permutation([2, 1, 4, 3])
    τ = Permutation([1, 3, 4, 2])
    x = 2
	S = [σ, τ]
	Δ, schreierVector = transversalSchreier(S, x)
	for (idx, δ) in pairs(Δ)
		δ == x && continue
		k = δ^inv(S[schreierVector[δ]])
		@test findfirst(==(k), Δ) < idx
	end
end

@testset "SchreierTransversal_2" begin
    σ = Permutation([2, 1, 4, 3])
    τ = Permutation([1, 3, 4, 2])
    x = 2
	S = [σ, τ]
	Δ, schreierVector = transversalSchreier(S, x)
	@test length(Δ) == 4
	for δ in Δ
		@test x^representative(δ, S, Δ, schreierVector) == δ
	end
end

