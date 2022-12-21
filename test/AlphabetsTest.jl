using ComputationalGroupTheory
using Test

@testset "Alphabets_1" begin
    A = Alphabet("a", "b", "c")
    @test length(A) == 3
    
    @test A["a"] == 1
    @test A["b"] == 2
    @test A["c"] == 3

    @test A[1] == "a"
    @test A[2] == "b"
    @test A[3] == "c"
end

@testset "Alphabets_2" begin
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