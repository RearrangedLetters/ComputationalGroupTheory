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

A = Σ"xXy"
setinverse!(A, "x", "X")

B = Σ"xXyYzZ"
setinverse!(B, "x", "X")
setinverse!(B, "y", "Y")
setinverse!(B, "z", "Z")

@testset "WordMacro_1" begin
    @test one(stringword"") == one(Word(["x"]))
    @test stringword"a" == Word(["a"])
    @test stringword"aAyx" == Word(["a", "A", "y", "x"])
end

@testset "FreeRewriteBitvector_1" begin
    @test freeRewriteBV!(stringword"", A) == stringword""

    @test freeRewriteBV!(stringword"xXxxX", A) == stringword"x"
    @test freeRewriteBV!(stringword"xxxxXXxXX", A) == stringword"x"
    @test freeRewriteBV!(stringword"xXxxX" * stringword"xxxxXXxXX", A) == stringword"xx"

    @test freeRewriteBV!(stringword"xxyzZYXX", B) == stringword""
    @test freeRewriteBV!(stringword"xaX", B) == stringword"xaX"
end

@testset "FreeRewrite" begin  # Lecture version & test
	A = Σ"xyX"  # CGT.Alphabet([:x, :y, :X])
	setinverse!(A, "x", "X")
	x, X, y = stringword"x", stringword"X", stringword"y"  # CGT.Word([A[:x]]), CGT.Word([A[:X]]), CGT.Word([A[:y]])
	@test rewrite(x*X, A) == one(x)
	@test rewrite(y*x*X, A) == y
	@test rewrite(X*y*x, A) == X*y*x

	@test rewrite(X*x*X, A) == X

	setinverse!(A, "y", "y")
	@test isone(rewrite(y*x*X*y, A))
end

@testset "Rewrite Algorithms Performance Comparison" begin
    letters = String[]
    for i ∈ 1:42000
        push!(letters, "x")
        push!(letters, "X")
    end
    w = Word(letters)
    #benchmark_BV = @benchmark freeRewriteBV!($w, A)
    println("buffer rewrite: ", @time freeRewriteBV!(w, A))
    w = Word(letters)
    #benchmark_lecture = @benchmark rewrite($w, A)
    println("lecture rewrite: ", @time rewrite(w, A))
    #println(judge(median(benchmark_BV), median(benchmark_lecture)))
end