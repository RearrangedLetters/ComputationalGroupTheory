using ComputationalGroupTheory
using Test
include("SmallPermutationGroups.jl")


@testset "OrderTests" begin
    @test order(I)  == 1
    @test order(C₂) == 2
    @test order(C₃) == 3
    @test order(C₄) == 4
    @test order(V)  == 4
    @test order(C₅) == 5
    @test order(C₆) == 6
    @test order(S₃) == 6
    @test order(C₇) == 7

    @test order(A₄) == 12
end

@testset "AbelianTest" begin
    @test isAbelian(I)  == true
    @test isAbelian(C₂) == true
    @test isAbelian(C₃) == true
    @test isAbelian(C₄) == true
    @test isAbelian(V)  == true
    @test isAbelian(C₅) == true
    @test isAbelian(C₆) == true
    @test isAbelian(S₃) == false
    @test isAbelian(C₇) == true

    @test isAbelian(A₄) == false
end