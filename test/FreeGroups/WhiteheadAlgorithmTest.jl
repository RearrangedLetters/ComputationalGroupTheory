using ComputationalGroupTheory
using Test
using ProfileView
using Random

Random.seed!(42)

const 𝟙    = Word(:𝟙)
const 𝟙⁻   = Word(:𝟙⁻)
const 𝟙𝟙   = Word(:𝟙, :𝟙)
const 𝟙𝟙⁻  = Word(:𝟙, :𝟙⁻)
const 𝟙⁻𝟙  = Word(:𝟙⁻, :𝟙)
const 𝟙⁻𝟙⁻ = Word(:𝟙⁻, :𝟙⁻)

const Y = Basis(symmetric_alphabet"ab")

const T = Alphabet(:𝟙, :𝟙⁻)
setinverse!(T, :𝟙, :𝟙⁻)
const X = Basis(T)
const idₓ = FreeGroupAutomorphism(X)

const Z = Basis(symmetric_alphabet"abcdefg")

const N = NielsenAutomorphisms(X)

@testset "Primitive elements in ℤ" begin
    #=
    Now F₁ ≅ ⟨𝟙⟩ ≅ ℤ, and there are exactly two primitive elements,
    namely 𝟙 ≙ 1 and -1 ≙ 𝟙⁻.
    =#
    @test isprimitive_nielsenfirst(𝟙, X)
    @test isprimitive_naive(𝟙, X)
    @test isprimitive_nielsenfirst(𝟙⁻, X)
    @test isprimitive_naive(𝟙⁻, X)

    found_automorphism₁, σ₁₁, σ₁₂, τ₁ = whitehead_nielsenfirst(𝟙, 𝟙⁻, X)
    @test found_automorphism₁
    @test length(σ₁₁) == length(σ₁₂) == 0
    
    found_automorphism₂, σ₂₁, σ₂₂, τ₂ = whitehead_naive(𝟙, 𝟙⁻, X)
    @test found_automorphism₂
    @test length(σ₂₁) == length(σ₂₂) == 0
    
    τ₁ = compose(τ₁)
    τ₂ = compose(τ₂)
    @test τ₁ == τ₂ == FreeGroupAutomorphism(X, [𝟙⁻])

    #=
    Now we assert that neither 2 ≙ 𝟙𝟙 nor -2 ≙ 𝟙⁻𝟙⁻ are primitive.
    =#
    @test !isprimitive_nielsenfirst(𝟙𝟙, X)
    @test !isprimitive_naive(𝟙𝟙, X)
    @test !isprimitive_nielsenfirst(𝟙⁻𝟙⁻, X)
    @test !isprimitive_naive(𝟙⁻𝟙⁻, X)
end

@testset "Primitivity (Naive)" begin
    @test !isprimitive_naive(word"", X)
    @test !isprimitive_naive(word"", Y)
    @test !isprimitive_naive(word"", Z)
    
    @test isprimitive_naive(word"a", Y)
    @test isprimitive_naive(word"ab", Y)
    
    @test !isprimitive_naive(word"ABab", Y)
    @test !isprimitive_naive(word"aa", Y)
    @test !isprimitive_naive(word"abaB", Y)
    
    @test isprimitive_naive(word"abcde", Z)
    @test isprimitive_naive(word"ABCD", Z)
end

@testset "Primitivity (Nielsen First)" begin
    @test !isprimitive_nielsenfirst(word"", X)
    @test !isprimitive_nielsenfirst(word"", Y)
    @test !isprimitive_nielsenfirst(word"", Z)
    
    @test isprimitive_nielsenfirst(word"a", Y)
    @test isprimitive_nielsenfirst(word"ab", Y)
    
    @test !isprimitive_nielsenfirst(word"ABab", Y)
    @test !isprimitive_nielsenfirst(word"aa", Y)
    @test !isprimitive_nielsenfirst(word"abaB", Y)

    @test isprimitive_nielsenfirst(word"abcde", Z)
    @test isprimitive_nielsenfirst(word"ABCDE", Z)
end

@testset "Primitivity (Nielsen Only)" begin
    @test !isprimitive_nielsenonly(word"", X)
    @test !isprimitive_nielsenonly(word"", Y)
    @test !isprimitive_nielsenonly(word"", Z)
    
    @test isprimitive_nielsenonly(word"a", Y)
    @test isprimitive_nielsenonly(word"ab", Y)
    
    @test !isprimitive_nielsenonly(word"ABab", Y)
    @test !isprimitive_nielsenonly(word"aa", Y)
    @test !isprimitive_nielsenonly(word"abaB", Y)

    @test isprimitive_nielsenonly(word"abcdefg", Z)
    @test isprimitive_nielsenonly(word"ABCDE", Z)
    @test isprimitive_nielsenonly(word"aaaabbbccd", Z)
end

@testset "Nielsen Graph (1)" begin
    G = AutomorphismGraph(X, wordlengths=[1], automorphisms=N, usecyclicwords=false)

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
    G = AutomorphismGraph(X, wordlengths=[2], automorphisms=N, usecyclicwords=false)

    @test order(G) == 2
    @test size(G)  == 2
    @test 𝟙 * 𝟙   ∈ G
    @test 𝟙 * 𝟙⁻  ∉ G
    @test 𝟙⁻ * 𝟙  ∉ G
    @test 𝟙⁻ * 𝟙⁻ ∈ G

    # There are no identities:
    @test length(edges(G, 𝟙𝟙,   𝟙𝟙))   == 0
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙⁻𝟙⁻)) == 0

    # The only edges are induced by the inversion, the only automorphism in N
    @test length(edges(G, 𝟙⁻𝟙⁻, 𝟙𝟙))  == 1
    @test length(edges(G, 𝟙𝟙, 𝟙⁻𝟙⁻))  == 1
end

@testset "Nielsen Graph (3)" begin
    M = NielsenAutomorphisms(Y)
    G = AutomorphismGraph(Y, wordlengths=[2], automorphisms=M, usecyclicwords=false)

    @test order(G) == 12
    @test size(G)  == 4 * 1 + 8 * 2
end

@testset "Nielsen Graph Path Test (1)" begin
    G = AutomorphismGraph(X, wordlengths=[1], automorphisms=N, usecyclicwords=false)

    τ₁ = connect_depthfirst(G, 𝟙, 𝟙)
    @test length(τ₁) == 0
    
    τ₂ = connect_depthfirst(G, 𝟙, 𝟙⁻)
    @test length(τ₂) == 1
    @test first(τ₂) == FreeGroupAutomorphism(X, [𝟙⁻])
    
    τ₃ = connect_depthfirst(G, 𝟙⁻, 𝟙)
    @test length(τ₃) == 1
    @test first(τ₃) == FreeGroupAutomorphism(X, [𝟙⁻])
    
    τ₄ = connect_depthfirst(G, 𝟙⁻, 𝟙⁻)
    @test length(τ₄) == 0
end

@testset "Nielsen Graph Path Test (2)" begin
    G = AutomorphismGraph(X, wordlengths=[3], automorphisms=N, usecyclicwords=false)

    τ₁ = connect_depthfirst(G, Word(:𝟙, :𝟙, :𝟙), Word(:𝟙⁻, :𝟙⁻, :𝟙⁻))
    @test length(τ₁) > 0
    @test compose(τ₁)(Word(:𝟙, :𝟙, :𝟙)) == Word(:𝟙⁻, :𝟙⁻, :𝟙⁻)
end

@testset "Nielsen Graph Path Test (2)" begin
    M = NielsenAutomorphisms(Y)
    G = AutomorphismGraph(Y, wordlengths=[2], automorphisms=M, usecyclicwords=false)

    τ₁ = connect_depthfirst(G, word"ab", word"BA")
    @test length(τ₁) == 0

    τ₂ = connect_depthfirst(G, word"bA", word"aa")
    @test length(τ₂) == 0
end

@testset "Automorphism Graph (1)" begin
    G = AutomorphismGraph(X, wordlengths=[1], usecyclicwords=false)
    # The vertices are 1 and -1
    @test order(G) == 2
    # The edges correspond to the inversions
    @test size(G)  == 2
end

@testset "Whitehead Word Reduction (1)" begin
    Y₂ = Basis(symmetric_alphabet"a")
    v₁ = word"a"
    for i ∈ 1:3
        w, σ, has_been_reduced = whitehead_reduce(v₁^i, Y₂)
        @test w == v₁^i
        @test isnothing(σ)
        @test !has_been_reduced
    end
end

@testset "Whitehead Word Reduction (2)" begin
    Y₂ = Basis(symmetric_alphabet"ab")

    v₁ = word"ab"
    w₁, σ₁, has_been_reduced₁ = whitehead_reduce(v₁, Y₂)
    @test w₁ == word"b"
    @test σ₁ == FreeGroupAutomorphism(Y₂, [word"a", word"Ab"])
    @test has_been_reduced₁
    
    v₂ = word"ba"
    w₂, σ₂, has_been_reduced₂ = whitehead_reduce(v₂, Y₂)
    @test w₂ == word"a"
    @test σ₂ == FreeGroupAutomorphism(Y₂, [word"Ba", word"b"])
    @test has_been_reduced₁
end

@testset "Inverse through Whitehead's Algorithm (1)" begin
    B₃  = Basis(symmetric_alphabet"abc")
    σ   = FreeGroupAutomorphism(B₃, [word"a", word"Aba", word"Aca"])
    σ⁻¹ = inv(σ)
    id₁ = compose(σ, σ⁻¹)
    id₂ = compose(σ⁻¹, σ)

    for x ∈ alphabet(B₃)
        @test id₁(x) == x
        @test id₂(x) == x
    end

    w = cyclically_reduce!(word"bAacBBbCbcAbCaAbaBBCAb", alphabet(B₃))
    @test id₁(w) == w
    @test id₂(w) == w
end

@testset "Compose (1)" begin
    τ₁ = FreeGroupAutomorphism(Y, [word"a", word"aab"])
    τ₂ = FreeGroupAutomorphism(Y, [word"aa", word"b"])
    @test compose([τ₁, τ₂]) == FreeGroupAutomorphism(Y, [word"aa", word"aab"])
    @test compose([τ₂, τ₁]) == FreeGroupAutomorphism(Y, [word"aa", word"aaaab"])
end

begin
    const X₃ = Basis(symmetric_alphabet"abc")
    
    const randomletters₁ = rand(collect('a':'c') ∪ collect('A':'C'), 10)
    const w₁ = cyclically_reduce(Word([Symbol(s) for s ∈ randomletters₁]), alphabet(X₃))
    
    const randomletters₂ = rand(collect('a':'c') ∪ collect('A':'C'), 10)
    const w₂ = cyclically_reduce(Word([Symbol(s) for s ∈ randomletters₂]), alphabet(X₃))

    const w₃ = Word([:b, :c, :A, :b, :a, :b, :A, :b])

    @testset "Minimize (1)" begin
        w, τ, isshorter = minimize!(deepcopy(w₁), X₃)
        @test length(τ) > 1
        @test apply(τ, w₁) == w
        @test isshorter
        @test length(w) < length(w₁)
    end

    @testset "Whitehead Naive (1)" begin
        success, σ₁, σ₂, τ = whitehead_naive(w₁, w₂, X₃)
        w₁′ = apply(σ₁, w₁)
        w₂′ = apply(σ₂, w₂)
        println(w₁′)
        println(w₂′)

        @test !success
        @test length(w₁′) != length(w₂′)
        @test !arecyclicallyequal(w₁′, w₂′)
    end

    @testset "Whitehead Nielsen First Heuristic (1)" begin
        success, σ₁, σ₂, τ = whitehead_nielsenfirst(w₁, w₂, X₃)
        w₁′ = apply(σ₁, w₁)
        w₂′ = apply(σ₂, w₂)

        @test !success
        @test length(w₁′) != length(w₂′)
        @test !arecyclicallyequal(w₁′, w₂′)
    end

    @testset "Whitehead Nielsen First Heuristic (2)" begin
        success, σ₁, σ₂, τ = whitehead_nielsenfirst(w₁, w₃, X₃)
        w₁′ = apply(σ₁, w₁)
        w₃′ = apply(σ₂, w₃)

        @test success
        @test length(w₁′) == length(w₃′)
        @test arecyclicallyequal(apply(τ, w₁′), w₃′)
    end
end