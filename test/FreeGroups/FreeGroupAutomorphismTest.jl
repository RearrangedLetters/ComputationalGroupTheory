using ComputationalGroupTheory
using Test

@testset "Whitehead Automorphisms (1)" begin
    X = Basis(symmetric_alphabet"a")
    W = WhiteheadAutomorphisms(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == 2
end

@testset "Whitehead Automorphisms (2)" begin
    X = Basis(symmetric_alphabet"ab")
    W = WhiteheadAutomorphisms(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 20
end

@testset "Whitehead Automorphisms Type I (1)" begin
    X = Basis(symmetric_alphabet"a")
    W = WhiteheadAutomorphismsTypeI(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == 2
    
    @test FreeGroupAutomorphism(X, [word"a"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"A"]) ∈ automorphisms
end

@testset "Whitehead Automorphisms Type I (2)" begin
    X = Basis(symmetric_alphabet"ab")
    W = WhiteheadAutomorphismsTypeI(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 8

    @test FreeGroupAutomorphism(X, [word"b", word"a"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"B", word"a"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"b", word"A"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"B", word"A"]) ∈ automorphisms

    @test FreeGroupAutomorphism(X, [word"a", word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"A", word"b"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"a", word"B"]) ∈ automorphisms
    @test FreeGroupAutomorphism(X, [word"A", word"B"]) ∈ automorphisms
end

@testset "Whitehead Automorphisms Type I (3)" begin
    X = Basis(symmetric_alphabet"abcde")
    W = WhiteheadAutomorphismsTypeI(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == factorial(length(X)) * 2^length(X) == length(W)
end

@testset "Whitehead Automorphisms Type II (1)" begin
    X = Basis(symmetric_alphabet"a")
    W = WhiteheadAutomorphismsTypeII(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 0
end

@testset "Whitehead Automorphisms Type II (2)" begin
    X = Basis(symmetric_alphabet"ab")
    W = WhiteheadAutomorphismsTypeII(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 12
end

@testset "Whitehead Automorphisms Type II (3)" begin
    X = Basis(symmetric_alphabet"abc")
    W = WhiteheadAutomorphismsTypeII(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 90
end

@testset "Whitehead Automorphisms Type II (4)" begin
    X = Basis(symmetric_alphabet"abcd")
    W = WhiteheadAutomorphismsTypeII(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W) == 504
end

@testset "Whitehead Automorphisms Type II (5)" begin
    X = Basis(symmetric_alphabet"abcde")
    W = WhiteheadAutomorphismsTypeII(X)
    automorphisms = Vector{FreeGroupAutomorphism}()
    for σ ∈ W push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == length(W)
end

@testset "FreeGroupAutomorphism Test" begin
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
    @test σ(word"aBbA") == ε
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
    A = Alphabet(:𝟙, :𝟙⁻)
    setinverse!(A, :𝟙, :𝟙⁻)
    X = Basis(A)

    automorphisms = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ NielsenAutomorphisms(X) push!(automorphisms, σ) end

    @test length(automorphisms) == length(Set(automorphisms))
    @test length(automorphisms) == 1
    σ = first(automorphisms)
    @test σ(:𝟙)  == :𝟙⁻
    @test σ(:𝟙⁻) == :𝟙
end

@testset "Verify NielsenAutomorphisms (2)" begin
    X = symmetric_alphabet"ab"
    N = NielsenAutomorphisms(Basis(X))

    automorphisms = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N
        push!(automorphisms, σ)
    end
    
    @test length(automorphisms) == length(Set(automorphisms))
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
    X₁ = Basis(symmetric_alphabet"abc")
    X₂ = Basis(symmetric_alphabet"abcdefg")
    N₁ = NielsenAutomorphisms(X₁)
    N₂ = NielsenAutomorphisms(X₂)

    automorphisms₁ = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N₁
        push!(automorphisms₁, σ)
    end
    @test length(automorphisms₁) == length(Set(automorphisms₁))
    @test length(automorphisms₁) == length(N₁)

    automorphisms₂ = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ N₂
        push!(automorphisms₂, σ)
    end
    @test length(automorphisms₂) == length(Set(automorphisms₂))
    @test length(automorphisms₂) == length(N₂)
end