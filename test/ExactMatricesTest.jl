using ComputationalGroupTheory
using Test

@testset "Echelonize_1" begin
    A = [1//1 0//1;
         0//1 1//1]
    B = echelonize(A)
    @test B == A
end

@testset "Echelonize_2" begin
    A = [2//1 0//1;
         0//1 2//1]
    B = echelonize(A)
    @test B == one(A)
end

@testset "Echelonize_3" begin
    A = [2//1 4//1;
         4//1 8//1]
    B = echelonize(A)
    @test B == [1//1 2//1;
                0//1 0//1]
end

@testset "Nullspace_1" begin
    A = [1 2 3; 2 4 6]
    A = convert(Matrix{Rational{BigInt}}, A)
    nullspace = ComputationalGroupTheory.nullspace(A)
    @test transpose(nullspace[1]) * A == transpose(zero(A[1, :]))
end

@testset "Nullspace_2" begin
    for i in 1:42
        A = rand(1:2//1, 10, 10)
        nullspace = ComputationalGroupTheory.nullspace(A)
        for k in nullspace
            @test transpose(k) * A == transpose(zero(A[1, :]))
        end
    end
end