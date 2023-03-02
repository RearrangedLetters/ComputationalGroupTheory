using ComputationalGroupTheory
using Test

#= @testset "Whitehead Automorphisms Type I (1)" begin
    X = Basis(symmetric_alphabet"a")
    W = WhiteheadAutomorphismsTypeI(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == 1
    
    σ = first(automorphisms)
    @test σ == FreeGroupAutomorphism(X, [word"a"])
end =#

@testset "Whitehead Automorphisms Type I (2)" begin
    X = Basis(symmetric_alphabet"ab")
    W = WhiteheadAutomorphismsTypeI(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == 2

    σ₁ = first(automorphisms)
    σ₂ = last(automorphisms)
    @test σ₁ == FreeGroupAutomorphism(X, [word"a", word"b"])
    @test σ₂ == FreeGroupAutomorphism(X, [word"b", word"a"])
end

@testset "Whitehead Automorphisms Type I (3)" begin
    X = Basis(symmetric_alphabet"abcde")
    W = WhiteheadAutomorphismsTypeI(X)
    count = 0
    for _ ∈ W count += 1 end

    @test count == factorial(length(X))
end

@testset "Whitehead Automorphisms Type II (1)" begin
    
end

@testset "Whitehead Automorphisms Type II (2)" begin
    
end

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
end

@testset "Verify NielsenAutomorphisms (1)" begin
    X = Alphabet(:𝟙, :𝟙⁻)
    setinverse!(X, :𝟙, :𝟙⁻)

    automorphisms = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ NielsenAutomorphisms(X) push!(automorphisms, σ) end
    @test length(automorphisms) == 1
    σ = first(automorphisms)
    @test σ(:𝟙)  == :𝟙⁻
    @test σ(:𝟙⁻) == :𝟙
end

@testset "Verify NielsenAutomorphisms (2)" begin
    X = symmetric_alphabet"ab"
    N = NielsenAutomorphisms(X)

    automorphisms = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N
        push!(automorphisms, σ)
    end
    
    @test length(automorphisms) == length(N)
    @test FreeGroupAutomorphism(X, [word"A",  word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"ba", word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"Ba", word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"ab", word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"aB", word"b"]) ∈ automorphisms
    
    @test FreeGroupAutomorphism(X, [word"a", word"B"])  ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"a", word"ab"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"a", word"Ab"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"a", word"ba"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"a", word"bA"]) ∈ automorphisms
end

@testset "Verify NielsenAutomorphisms (3)" begin
    X₁ = symmetric_alphabet"abc"
    X₂ = symmetric_alphabet"abcdefg"
    N₁ = NielsenAutomorphisms(X₁)
    N₂ = NielsenAutomorphisms(X₂)

    automorphisms₁ = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N₁
        push!(automorphisms₁, σ)
    end
    @test length(automorphisms₁) == length(N₁)

    automorphisms₂ = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N₂
        push!(automorphisms₂, σ)
    end
    @test length(automorphisms₂) == length(N₂)
end =#

#= @testset "Count Free Group Automorphisms 1" begin
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
end =#

#= @testset "Count Free Group Automorphisms (2)" begin
    X = symmetric_alphabet"a"
    @info length(WhiteheadAutomorphisms(X))
    W = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ WhiteheadAutomorphisms(X)
        push!(W, σ)
        @info σ
    end
    
end =#

