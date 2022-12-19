using ComputationalGroupTheory
using Test

# Based on en.wikipedia.org/wiki/List_of_small_groups

I = [perm"()"]
S₁ = I
A₂ = I

C₂ = [perm"(2, 1)"]
S₂ = C₂
D₂ = C₂

C₃ = [perm"(2, 3, 1)"]
A₃ = C₃

C₄ = [perm"(2, 3, 4, 1)"]
V = [perm"(1, 2)(3, 4)", perm"(1, 3)(2, 4)"]
K₄ = V

C₅ = [perm"(2, 3, 4, 5, 1)"]

C₆ = [perm"(2, 3, 4, 5, 6, 1)"]
S₃ = [perm"(1, 2)", perm"(2, 3)"]
D₆ = S₃

C₇ = [perm"(2, 3, 4, 5, 6, 7, 1)"]


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
end