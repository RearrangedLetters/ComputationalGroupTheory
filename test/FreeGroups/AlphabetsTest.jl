using ComputationalGroupTheory
using Test

@testset "Alphabets 1" begin
    A = Alphabet("a", "b", "c")
    @test length(A) == 3
    
    @test A["a"] == 1
    @test A["b"] == 2
    @test A["c"] == 3

    @test A[1] == "a"
    @test A[2] == "b"
    @test A[3] == "c"
end

@testset "Alphabets 2" begin
    A = Alphabet(["x", "y", "z"])
    for a ∈ A 
        @test !hasinverse(A, a)
    end
    setinverse!(A, "x", "X")
    setinverse!(A, "y", "Y")
    setinverse!(A, "z", "Z")
    for a ∈ A
        @test hasinverse(A, a)
    end
end

@testset "Alphabets 3" begin
    A = symmetric_alphabet"abc"
    @test isinverse(A, :a, :A)
    @test isinverse(A, :B, :b)
    @test !isinverse(A, :a, :B)
    @test !isinverse(A, :c, :c)
end