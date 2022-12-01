using ComputationalGroupTheory
using Test

@testset "Backtrack_1" begin
    a = perm"(1, 2, 3, 4)"
    b = perm"(2, 3)"
    S = [a, b]
    @info enumerateGroup(S)
    @info length(enumerateGroup(S))
end