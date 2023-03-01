using ComputationalGroupTheory
using Test

const X = Alphabet(:𝟙, :𝟙⁻)
setinverse!(X, :𝟙, :𝟙⁻)
const 𝟙    = Word(:𝟙)
const 𝟙⁻   = Word(:𝟙⁻)
const 𝟙𝟙   = Word(:𝟙, :𝟙)
const 𝟙𝟙⁻  = Word(:𝟙, :𝟙⁻)
const 𝟙⁻𝟙  = Word(:𝟙⁻, :𝟙)
const 𝟙⁻𝟙⁻ = Word(:𝟙⁻, :𝟙⁻)

const Y = symmetric_alphabet"ab"
const a = word"a"
const A = word"A"
const b = word"b"
const B = word"B"


#= @testset "Nielsen Graph (1)" begin
    N = NielsenAutomorphisms(X)
    G = AutomorphismGraph(X, wordlength=1, automorphisms=N)

    @test order(G) == 2
    @test size(G)  == 2
    @test 𝟙  ∈ G
    @test 𝟙⁻ ∈ G

    @test length(edges(G, 𝟙, 𝟙))   == 0
    @test length(edges(G, 𝟙, 𝟙⁻))  == 1
    @test length(edges(G, 𝟙⁻, 𝟙))  == 1
    @test length(edges(G, 𝟙⁻, 𝟙⁻)) == 0
end

@testset "Nielsen Graph (2)" begin
    N = NielsenAutomorphisms(X)
    G = AutomorphismGraph(X, wordlength=2, automorphisms=N)

    @test order(G) == 4
    @test size(G)  == 4
    @test 𝟙 * 𝟙   ∈ G
    @test 𝟙 * 𝟙⁻  ∈ G
    @test 𝟙⁻ * 𝟙  ∈ G
    @test 𝟙⁻ * 𝟙⁻ ∈ G

    # There are no identities:
    @test length(edges(G, 𝟙𝟙,   𝟙𝟙))   == 0
    @test length(edges(G, 𝟙⁻𝟙,  𝟙⁻𝟙))  == 0
    @test length(edges(G, 𝟙𝟙⁻,  𝟙𝟙⁻))  == 0
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙⁻𝟙⁻)) == 0

    # The only edges are induced by the inversion, the only automorphism in N
    @test length(edges(G, 𝟙𝟙, 𝟙⁻𝟙⁻))  == 1
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙𝟙))  == 1
    @test length(edges(G, 𝟙⁻𝟙, 𝟙𝟙⁻))  == 1
    @test length(edges(G, 𝟙𝟙⁻, 𝟙⁻𝟙))  == 1
    
    # The remaining possibilities:
    @test length(edges(G, 𝟙𝟙, 𝟙⁻𝟙))   == 0
    @test length(edges(G, 𝟙𝟙, 𝟙𝟙⁻))   == 0
    @test length(edges(G, 𝟙⁻𝟙, 𝟙𝟙))   == 0
    @test length(edges(G, 𝟙⁻𝟙, 𝟙⁻𝟙⁻)) == 0
    @test length(edges(G, 𝟙𝟙⁻, 𝟙𝟙))   == 0
    @test length(edges(G, 𝟙𝟙⁻, 𝟙⁻𝟙⁻)) == 0
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙𝟙⁻)) == 0
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙⁻𝟙)) == 0
end =#

@testset "Nielsen Graph (3)" begin
    N = NielsenAutomorphisms(Y)
    G = AutomorphismGraph(Y, wordlength=2, automorphisms=N)

    @test order(G) == 4^2
    @test size(G)  == order(G) * length(N)
end


@testset "Automorphism Graph (1)" begin
    G = AutomorphismGraph(X, wordlength=1)
    @info G
    # The vertices are 1 and -1
    @test order(G) == 2
    # The edges are the two identity self-loops and the inversion
    @test size(G)  == 3
end

@testset "Whitehead Word Reduction (1)" begin
    Y = symmetric_alphabet"a"
    v₁ = word"a"
    for i ∈ 1:3
        w, σ, has_been_reduced = whitehead_reduce!(Y, v₁^i)
        @test w == v₁
        @test isnothing(σ)
        @test !has_been_reduced
    end
end

@testset "Whitehead Word Reduction (2)" begin
    Y = symmetric_alphabet"ab"

    v₁ = word"ab"  # can be reduced by setting x = a, y = b and x ↦ aB
    w₁, σ₁, has_been_reduced₁ = whitehead_reduce!(Y, v₁)
    @test w₁ == word"a"
    @test σ₁ == FreeGroupAutomorphism(Y, [word"aB", word"b", word"bA", word"B"])
    @test has_been_reduced₁
    
    v₂ = word"ba"  # can be reduced by setting x = b, y = a and x ↦ bA
    w₂, σ₂, has_been_reduced₂ = whitehead_reduce!(Y, v₂)
    @test w₂ == word"a"
    @test σ₂ == FreeGroupAutomorphism(Y, [word"a", word"bA", word"A", word"aB"])
    @test has_been_reduced₁
end

@testset "Automorphism Graph (2)" begin
    Y = symmetric_alphabet"ab"

    G = AutomorphismGraph(Y, wordlength=2)
    # Vertices:
    #   :a:a, :a:A, :a:b, :a:B,
    #   :b:a, :b:A, :b:b, :b:B
    #   :A:a, :A:A, :A:b, :A:B,
    #   :B:a, :B:A, :B:b, :B:B,
    @test order(G) == 16

    @test size(G)  == 16 + k
end

@testset "Automorphism Graph (3)" begin
    Y = symmetric_alphabet"abc"

    G = AutomorphismGraph(Y, wordlength=3)
    @test order(G) == 1
    @test size(G)  == 1
end

#=
@testset "Primitive elements in ℤ" begin
    X = Alphabet(:𝟙)
    setinverse!(X, :𝟙, :𝟙⁻)

    @test freeRewriteBV!(word"𝟙𝟙⁻", X) == word""

    #=
    Now F₁ ≅ ⟨𝟙⟩ ≅ ℤ, and there are exactly two primitive elements,
    namely 𝟙 ≙ 1 and -1 ≙ 𝟙⁻.
    =#
    @test isprimitive_naive(X, word"𝟙")
    @test isprimitive_naive(X, word"𝟙⁻")

    #=
    Now we assert that neither 2 ≙ 𝟙𝟙 nor -2 ≙ 𝟙⁻𝟙⁻ are primitive.
    =#
    @test !isprimitive_naive(X, word"𝟙𝟙")
    @test !isprimitive_naive(X, word"𝟙⁻𝟙⁻")
end =#