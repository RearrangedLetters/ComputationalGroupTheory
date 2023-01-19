using ComputationalGroupTheory
using Test

@testset "LenLex" begin
	A = Alphabet([:a, :A, :b, :B])
    setinverse!(A, :a, :A)
    setinverse!(A, :b, :B)

    ord = LenLex(A, [:a, :A, :b, :B])

    @test ord isa Base.Order.Ordering

    u1 = Word([1,2])
    u3 = Word([1,3])
    u4 = Word([1,2,3])
    u5 = Word([1,4,2])

    @test lt(ord, u1, u3) == true
    @test lt(ord, u3, u1) == false
    @test lt(ord, u3, u4) == true
    @test lt(ord, u4, u5) == true
    @test lt(ord, u5, u4) == false
    @test lt(ord, u1, u1) == false
end