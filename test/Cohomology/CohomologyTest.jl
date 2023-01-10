using ComputationalGroupTheory
using Test

Λ₁ = [1  0  0 -1  0  0;
      0  1  0  0 -1  0;
      0  0  1  0  0 -1;
      0  0  0 -1  0  0;
      0  0  0  0 -1  0;
      0  0  0  0  0 -1]

Λ₂ = [1  0  0  0  0  0;
      0  0  0  0 -1  0;
      0  0 -1  0  0  1;
      0  0  0  1  0  0;
      0  1  0  0 -1  0;
      0  0 -1  0  0  0]

X = [1, 2]  # [f = 1, g = 2]
R = [[(1, 1), (1, 1), (1, 1)],
     [(2, 1), (2, 1)],
     [(2, 1), (1, 1), (2, 1), (1, 1)]]
Λ = [Matrix{Rational{Int}}(Λ₁), Matrix{Rational{Int}}(Λ₂)]

@testset "H¹_1" begin
    ZM = Z1(X, R, Λ)
    Z¹ = nullspace(ZM)
    BM = B1(X, R, Λ)
    B¹ = rowspace(BM)
    @info ZM
    @info Z¹
    @info BM
    @info B¹
end

#=
This is Example 7.4 from the "Handbook of Computational Group Theory" on p. 249.
=#

# G = ⟨x, y | x³, y³, (x⁻¹y)²⟩ ≅ A₄ acts on the permutation module of 𝔽₃,
# the action matrices are:
𝟘 = Residue(0, 3)
𝟙 = Residue(1, 3)
Λ₁ = [𝟘 𝟘 𝟙 𝟘;
      𝟙 𝟘 𝟘 𝟘;
      𝟘 𝟙 𝟘 𝟘;
      𝟘 𝟘 𝟙 𝟘]

Λ₂ = [𝟙 𝟘 𝟘 𝟘;
      𝟘 𝟘 𝟙 𝟘;
      𝟘 𝟘 𝟘 𝟙;
      𝟘 𝟙 𝟘 𝟘]

X = [1, 2]
R = [[(1, 1), (1, 1), (1, 1)],
     [(2, 1), (2, 1), (2, 1)],
     [(1, -1), (2, 1), (1, -1), (2, 1)]]
    
Λ = [Λ₁, Λ₂]

r = length(X)
s = length(R)
d = size(Λ[1], 2)
zero1 = [[Residue(0, 3) for _ in 1:(d * r)] for _ in 1:(d * s)]
zero2 = [[Residue(0, 3) for _ in 1:(d * r)] for _ in 1:(d * s)]


@testset "H¹_2" begin
    ZM = Z1!(X, R, Λ, zero1)
    Z¹ = nullspace(ZM)
    BM = B1!(X, R, Λ, zero1)
    B¹ = rowspace(BM)
    @info ZM
    @info Z¹
    @info BM
    @info B¹
end