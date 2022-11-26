using Test
using ComputationalGroupTheory

@testset "Transversal_1" begin
    σ = Permutation([1, 3, 4, 2])
    τ = Permutation([1, 2, 4, 5, 3])
    x = 2
	T = Transversal([σ, τ], x)
	for δ in T.Ωᴳ
		@test x^T[δ] == δ
	end
end

@testset "FactoredTransversal_1" begin
    σ = Permutation([1, 3, 4, 2])
    τ = Permutation([1, 2, 4, 5, 3])
    x = 2
	T = FactoredTransversal([σ, τ], x)
	for δ in T.Ωᴳ
		@test x^T[δ] == δ
	end
end