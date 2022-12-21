using Test

@testset "PermutationGroup_1" begin
    G = PermutationGroup(perm"(1, 2, 3)", perm"(2, 3)")

    @test one(G) ∈ G
    @test one(G) isa eltpe(G)

    @test perm"(1, 2, 3)" ∈ G
    @test order(G) == 3  # todo: check order by hand, likely is not actually 3
    # check ∉
    # check for n random elements from the groups that they are in the group
    # check for multiples g * x if they are i nthe group
    # Random.seed!(42)
end

@test "PermutationGroupBasis_1" begin
    G = PermutationGroup(perm"(1, 2, 3)", perm"(2, 3)")

    g = 

    @test_throws permutationFromImages(stabilizerChain(G), )
end