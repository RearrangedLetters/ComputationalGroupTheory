using ComputationalGroupTheory
using Test

@testset "Residue_1" begin
    A = Residue(2, 5)
    B = Residue(3, 5)

    @test -A isa Residue
    @test A - B isa Residue
    @test zero(A) isa Residue
    @test A == A
    @test A != B
    @test A + zero(A) == A
    @test zero(A) - B == -B
    @test A - A == zero(A)
    @test -A + A == zero(A)
    @test A + 3 == zero(A)
    @test 3 + A == zero(A)
    @test big(3) + A == zero(A)
    @test A * B isa Residue
    @test A * inv(B) isa Residue
end