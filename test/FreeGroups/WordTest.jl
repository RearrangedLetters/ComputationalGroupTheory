using ComputationalGroupTheory
using Test
using BenchmarkTools

@testset "Test iterating over all words of fixed length" begin
    A = alphabet"ab"
    words₁ = Vector{Word{Symbol}}()
    for word ∈ Words(A, 1) push!(words₁, word) end
    @test words₁ == [word"a", word"b"]
    words₂ = Vector{Word{Symbol}}()
    for word ∈ Words(A, 2) push!(words₂, word) end
    @test words₂ == [word"aa", word"ba",
                     word"ab", word"bb"]
end

A = Alphabet(:x, :X, :y)
setinverse!(A, :x, :X)

B = symmetric_alphabet"xyz"

@testset "Word Macro" begin
    @test word"" == Word{Symbol}()
    @test word"a" == Word([:a])
    @test word"aAyx" == Word([:a, :A, :y, :x])
end

@testset "FreeRewrite with Bitvector" begin
    @test freerewriteBV!(word"", A) == word""

    @test freerewriteBV!(word"xXxxX", A) == word"x"
    @test freerewriteBV!(word"xxxxXXxXX", A) == word"x"
    @test freerewriteBV!(word"xXxxX" * word"xxxxXXxXX", A) == word"xx"

    @test freerewriteBV!(word"xxyzZYXX", B) == word""
    @test freerewriteBV!(word"xaX", B) == word"xaX"
end

@testset "FreeRewrite" begin  # Lecture version & test
	x, X, y = word"x", word"X", word"y"  # CGT.Word([A[:x]]), CGT.Word([A[:X]]), CGT.Word([A[:y]])
	@test rewrite(x*X, A) == one(x)
	@test rewrite(y*x*X, A) == y
	@test rewrite(X*y*x, A) == X*y*x

	@test rewrite(X*x*X, A) == X

	setinverse!(A, :y, :y)
	@test isone(rewrite(y*x*X*y, A))
end

@testset "Rewrite Algorithms Performance Comparison" begin
    letters = Symbol[]
    for i ∈ 1:420
        push!(letters, :x)
        push!(letters, :X)
    end
    w = Word(letters)
    #benchmark_BV = @benchmark freerewriteBV!($w, A)
    println("buffer rewrite: ")
    t₁ = @time freerewriteBV!(w, A)
    w = Word(letters)
    #benchmark_lecture = @benchmark rewrite($w, A)
    println("lecture rewrite: ")
    t₂ = @time rewrite(w, A)
    #println(judge(median(benchmark_BV), median(benchmark_lecture)))
end

@testset "Cyclic Reduction" begin
    A = symmetric_alphabet"ab"

    @test cyclically_reduce(word"aba", A) == word"aba"
    @test cyclically_reduce(word"Aba", A) == word"b"
    @test cyclically_reduce(word"bBabAAa", A) == word"b"
    @test cyclically_reduce(word"abBA", A) == word""
end