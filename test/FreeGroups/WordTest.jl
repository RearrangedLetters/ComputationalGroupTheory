using ComputationalGroupTheory
using Test
using BenchmarkTools

A = Σ"xXy"
setinverse!(A, "x", "X")

B = Σ"xXyYzZ"
setinverse!(B, "x", "X")
setinverse!(B, "y", "Y")
setinverse!(B, "z", "Z")

@testset "WordMacro_1" begin
    @test one(w"") == one(Word(["x"]))
    @test w"a" == Word(["a"])
    @test w"aAyx" == Word(["a", "A", "y", "x"])
end

@testset "FreeRewriteBitvector_1" begin
    @test freeRewriteBV!(w"", A) == w""

    @test freeRewriteBV!(w"xXxxX", A) == w"x"
    @test freeRewriteBV!(w"xxxxXXxXX", A) == w"x"
    @test freeRewriteBV!(w"xXxxX" * w"xxxxXXxXX", A) == w"xx"

    @test freeRewriteBV!(w"xxyzZYXX", B) == w""
    @test freeRewriteBV!(w"xaX", B) == w"xaX"
end

@testset "FreeRewrite" begin  # Lecture version & test
	A = Σ"xyX"  # CGT.Alphabet([:x, :y, :X])
	setinverse!(A, "x", "X")
	x, X, y = w"x", w"X", w"y"  # CGT.Word([A[:x]]), CGT.Word([A[:X]]), CGT.Word([A[:y]])
	@test rewrite(x*X, A) == one(x)
	@test rewrite(y*x*X, A) == y
	@test rewrite(X*y*x, A) == X*y*x

	@test rewrite(X*x*X, A) == X

	setinverse!(A, "y", "y")
	@test isone(rewrite(y*x*X*y, A))
end

@testset "Rewrite Algorithms Performance Comparison" begin
    letters = String[]
    for i ∈ 1:420000
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