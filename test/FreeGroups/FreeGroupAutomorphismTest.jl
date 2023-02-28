using ComputationalGroupTheory
using Test

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