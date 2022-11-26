using Test
using ComputationalGroupTheory

@testset "FactoredTransversal_1" begin
    σ = Permutation([1,3,4,2])
    τ = Permutation([1,2,4,5,3])
    x = 2
	Δ, T = transversalFactored([σ, τ], x)
	@test length(Δ) == 4
	for δ in Δ
		@test x^prod(T[δ]) == δ
	end
	Δ, T
end