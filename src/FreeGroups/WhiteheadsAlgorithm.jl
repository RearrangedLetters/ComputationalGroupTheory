include("Word.jl")
include("Alphabet.jl")
using Test
using Combinatorics
using Graphs

"""
The type Word can be used like a cyclic word by using getcyclicindex instead
of getindex. To avoid bugs that are hard to find, the usually indexing is not
cyclic by default.
"""
begin
    w = word"abba"
    cyclicword = Vector{Symbol}()
    for i ∈ 1:8
        push!(cyclicword, getcyclicindex(w, i))
    end
    @test cyclicword == word"abbaabba"
end

struct WhiteheadAutomorphisms
    rewritingSystem::RewritingSystem
    w::Word
    word_length::Int
    number_of_permutations::BigInt
    number_of_subsets::BigInt

    function WhiteheadAutomorphisms(rewritingSystem::RewritingSystem, w::Word)
        word_length = Base.length(w)
        new(rewritingSystem, w, word_length, factorial(big(word_length)), big(2)^word_length)
    end
end

function iterate(W::WhiteheadAutomorphisms)
    return Base.empty(W.rewritingSystem), (1, iterate(powerset(W.w.letters)))
end

function iterate(W::WhiteheadAutomorphisms, state)
    A = alphabet(W.rewritingSystem)
    i, powerset, j = state
    powerset_iterator = iterate(powerset, j)
    if i ≤ W.number_of_permutations
        σ = nthperm(letters(A), i)
        if !isnothing(powerset_iterator)
            subset_indices, _ = powerset_iterator
            for index ∈ subset_indices
                push!(rules, (Word(A[index]), Word(σ[index])))
            end
        else
            return iterate(W, (i + 1, powerset(letters(W.w.letters))))
        end
    else
        return nothing
    end
end