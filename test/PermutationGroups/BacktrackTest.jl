using ComputationalGroupTheory
using Test
include("SmallPermutationGroups.jl")

#= @testset "Backtrack_1" begin
    a = perm"(1, 2, 3, 4)"
    b = perm"(2, 3)"
    S = [a, b]
    @test length(backtrackRecursive(S)) == order(S)
end

@testset "DisplayBacktrackTree_1" begin
    FullA₄ = backtrackRecursive(A₄)
    @test length(FullA₄) == order(A₄)
    @test FullA₄ ⊆ Set(FullA₄)
end =#

@testset "IterativeBacktrack_1" begin
    # backtrack(A₄)
    for g ∈ A₄
        println(g)
    end
end