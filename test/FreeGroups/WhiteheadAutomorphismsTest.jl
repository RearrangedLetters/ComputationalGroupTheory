using ComputationalGroupTheory
using Test

#= @testset "Count Free Group Automorphisms 1" begin
#=
First we assert that the iteration protocol actually does the desired number
of iterations. This number is not equal to the actual number of Whitehead
automorphisms because the Nielsen automorphisms are covered twice.
=#
X₁ = symmetric_alphabet"a"
number_of_automorphisms = length(WhiteheadAutomorphisms(X₁))
@test number_of_automorphisms == 2
whitehead_automorphisms = collect(WhiteheadAutomorphisms(X₁))
id = whitehead_automorphisms[1]
σ  = whitehead_automorphisms[2]
@test id == FreeGroupAutomorphism(X₁, [word"a", word"A"])
@test σ  == FreeGroupAutomorphism(X₁, [word"A", word"a"])
end =#

@testset "Count Free Group Automorphisms (2)" begin
    X = symmetric_alphabet"a"
    @info length(WhiteheadAutomorphisms(X))
    W = Vector{FreeGroupAutomorphism{Symbol}}()
    for σ ∈ WhiteheadAutomorphisms(X)
        push!(W, σ)
        @info σ
    end
    
end

#= @testset "Verify Construction of Nielsen automorphisms" begin
    X₁ = symmetric_alphabet"a"
    σ = nielsen(X₁, 1, 1, 1)
    @test σ(:a) == word"A"
    @test σ(:A) == word"a"


    X₂ = symmetric_alphabet"ab"
    # In these tests x = a and y = b:
    @test nielsen(X₂, 1, 1, 1) == FreeGroupAutomorphism(X₂, [word"A",  word"b", word"a",  word"B"])
    @test nielsen(X₂, 1, 2, 2) == FreeGroupAutomorphism(X₂, [word"ba", word"b", word"AB", word"B"])
    @test nielsen(X₂, 1, 2, 3) == FreeGroupAutomorphism(X₂, [word"Ba", word"b", word"Ab", word"B"])
    @test nielsen(X₂, 1, 2, 4) == FreeGroupAutomorphism(X₂, [word"ab", word"b", word"BA", word"B"])
    @test nielsen(X₂, 1, 2, 5) == FreeGroupAutomorphism(X₂, [word"aB", word"b", word"bA", word"B"])
    @test isnothing(nielsen(X₂, 1, 2, 6))

    # In these tests x = b and y = a:
    @test nielsen(X₂, 2, 1, 2) == FreeGroupAutomorphism(X₂, [word"a", word"ab", word"A", word"BA"])
    @test nielsen(X₂, 2, 1, 3) == FreeGroupAutomorphism(X₂, [word"a", word"Ab", word"A", word"Ba"])
    @test nielsen(X₂, 2, 1, 4) == FreeGroupAutomorphism(X₂, [word"a", word"ba", word"A", word"AB"])
    @test nielsen(X₂, 2, 1, 5) == FreeGroupAutomorphism(X₂, [word"a", word"bA", word"A", word"aB"])
    @test isnothing(nielsen(X₂, 2, 1, 6))
end =#