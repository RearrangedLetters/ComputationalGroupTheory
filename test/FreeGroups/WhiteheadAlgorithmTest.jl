using ComputationalGroupTheory
using Test

#= @testset "FreeGroupAutomorphism Test" begin
    X = symmetric_alphabet"abc"
    ε = Word{Symbol}()

    id = FreeGroupAutomorphism(X)
    @test typeof(id(ε)) == Word{Symbol}
    @test typeof(id(:a)) == Word{Symbol}
    @test typeof(id(word"b")) == Word{Symbol}
    @test typeof(id(word"cc")) == Word{Symbol}

    @test id(ε) == ε
    @test id(:a) == :a
    @test id(:a) == word"a"
    @test id(word"a") == :a
    @test id(word"a") == word"a"
    @test id(word"c") == :c
    @test id(word"c") == word"c"
    @test id(word"abc") == word"abc"
    @test id(word"cab") == word"cab"

    @test id(:A) == :A
    @test id(:A) == word"A"
    @test id(word"A") == :A
    @test id(word"A") == word"A"
    @test id(word"C") == :C
    @test id(word"C") == word"C"
    @test id(word"ABC") == word"ABC"
    @test id(word"CAB") == word"CAB"

    σ = FreeGroupAutomorphism(X, [word"a", word"c", word"b"])
    @test σ(ε) == ε
    @test σ(:a) == :a
    @test σ(:b) == :c
    @test σ(:c) == :b
    @test σ(word"abc") == word"acb"
    @test σ(word"aaa") == word"aaa"
    @test σ(word"bb") == word"cc"
    @test σ(word"cccc") == word"bbbb"
    @test σ(word"bcbc") == word"cbcb"

    @test σ(:A) == :A
    @test σ(:B) == :C
    @test σ(:C) == :B
    @test σ(word"ABC") == word"ACB"
    @test σ(word"AAA") == word"AAA"
    @test σ(word"BB") == word"CC"
    @test σ(word"CCCC") == word"BBBB"
    @test σ(word"BCBC") == word"CBCB"
end =#

#= @testset "Verify Construction of Nielsen automorphisms" begin
    X₁ = symmetric_alphabet"a"
    σ = nielsen(X₁, 1, 1, 1)
    @test σ(:a) == word"A"
    @test σ(:A) == word"a"

    X₂ = symmetric_alphabet"ab"
    # In these tests x = a and y = b:
    @test nielsen(X₂, 1, 1, 1) == FreeGroupAutomorphism(X₂, [word"A",  word"b", word"a",  word"B"])
    @test nielsen(X₂, 1, 2, 2) == FreeGroupAutomorphism(X₂, [word"ba", word"b", word"AB", word"B"])
    @test nielsen(X₂, 1, 2, 3) == FreeGroupAutomorphism(X₂, [word"Ba", word"b", word"Ab", word"B"])
    @test nielsen(X₂, 1, 2, 4) == FreeGroupAutomorphism(X₂, [word"ab", word"b", word"BA", word"B"])
    @test nielsen(X₂, 1, 2, 5) == FreeGroupAutomorphism(X₂, [word"aB", word"b", word"bA", word"B"])
    @test isnothing(nielsen(X₂, 1, 2, 6))

    # In these tests x = b and y = a:
    @test nielsen(X₂, 2, 1, 2) == FreeGroupAutomorphism(X₂, [word"a", word"ab", word"A", word"BA"])
    @test nielsen(X₂, 2, 1, 3) == FreeGroupAutomorphism(X₂, [word"a", word"Ab", word"A", word"Ba"])
    @test nielsen(X₂, 2, 1, 4) == FreeGroupAutomorphism(X₂, [word"a", word"ba", word"A", word"AB"])
    @test nielsen(X₂, 2, 1, 5) == FreeGroupAutomorphism(X₂, [word"a", word"bA", word"A", word"aB"])
    @test isnothing(nielsen(X₂, 2, 1, 6))
end =#

@testset "Count Free Group Automorphisms (1)" begin
    #=
    First we assert that the iteration protocol actually does the desired number
    of iterations. This number is not equal to the actual number of Whitehead
    automorphisms because the Nielsen automorphisms are covered twice.
    =#
    X₁ = symmetric_alphabet"a"
    number_of_automorphisms = length(WhiteheadAutomorphisms(X₁))
    @test number_of_automorphisms == 2
    whitehead_automorphisms = collect(WhiteheadAutomorphisms(X₁))
    id = whitehead_automorphisms[1]
    σ  = whitehead_automorphisms[2]
    @test id == FreeGroupAutomorphism(X₁, [word"a", word"A"])
    @test σ  == FreeGroupAutomorphism(X₁, [word"A", word"a"])
end

@testset "Count Free Group Automorphisms (2)" begin
    X₂ = symmetric_alphabet"ab"
    println(length(WhiteheadAutomorphisms(X₂)))
end

@testset "Whitehead Word Reduction (1)" begin
    X = symmetric_alphabet"a"
    v₁ = word"a"
    for i ∈ 1:3
        w, σ, has_been_reduced = whitehead_reduce!(X, v₁^i)
        @test w == v₁
        @test isnothing(σ)
        @test !has_been_reduced
    end
end

@testset "Whitehead Word Reduction (2)" begin
    X = symmetric_alphabet"ab"

    v₁ = word"ab"  # can be reduced by setting x = a, y = b and x ↦ aB
    w₁, σ₁, has_been_reduced₁ = whitehead_reduce!(X, v₁)
    @test w₁ == word"a"
    @test σ₁ == FreeGroupAutomorphism(X, [word"aB", word"b", word"bA", word"B"])
    @test has_been_reduced₁
    
    v₂ = word"ba"  # can be reduced by setting x = b, y = a and x ↦ bA
    w₂, σ₂, has_been_reduced₂ = whitehead_reduce!(X, v₂)
    @test w₂ == word"a"
    @test σ₂ == FreeGroupAutomorphism(X, [word"a", word"bA", word"A", word"aB"])
    @test has_been_reduced₁
end

@testset "Automorphism Graph (1)" begin
    X = symmetric_alphabet"a"

    G = AutomorphismGraph(X, wordlength=1)
    # The vertices are :a and :A
    @test order(G) == 2
    # The edges are the two identity self-loops and the inversion
    @test size(G)  == 3
end

@testset "Automorphism Graph (2)" begin
    X = symmetric_alphabet"ab"

    G = AutomorphismGraph(X, wordlength=2)
    # Vertices:
    #   :a:a, :a:A, :a:b, :a:B,
    #   :b:a, :b:A, :b:b, :b:B
    #   :A:a, :A:A, :A:b, :A:B,
    #   :B:a, :B:A, :B:b, :B:B,
    @test order(G) == 16

    @test size(G)  == 16 + k
end

@testset "Automorphism Graph (3)" begin
    X = symmetric_alphabet"abc"

    G = AutomorphismGraph(X, wordlength=3)
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