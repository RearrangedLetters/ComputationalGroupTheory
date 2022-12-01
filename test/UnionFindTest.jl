using Test
using ComputationalGroupTheory
using Base.Iterators

@testset "UnionFind_1" begin
    n = 1000
    unionFind = UnionFind(n)
    @test length(collectBlocks(unionFind)) == n
    for i in 1:(n - 1)
        union(unionFind, i, i + 1)
    end
    blocks = collectBlocks(unionFind)
    @test length(blocks) == 1
    block = blocks[1]
    @test length(block) == n
    @test 1:n ⊆ block && block ⊆ 1:n
end

@testset "UnionFind_2" begin
    n = 666
    unionFind = UnionFind(n)
    for i in rand(1:n, convert(Int, floor(n / 4)))
        union!(unionFind, i, mod1(i * 2, n))
    end
    blocks = collectBlocks(unionFind)
    set = collect(Iterators.flatten(blocks))
    @test 1:n ⊆ set && set ⊆ 1:n
    @test Vector{UInt}() ∉ blocks
end