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

@testset "Count Whitehead Automorphisms" begin
    W₁ = WhiteheadAutomorphisms(rewritingSystem, word"a")
    length_W₁ = 0
    for _ in W₁ length_W₁ += 1 end
    @test length_W₁ == 2

    v = word"a"

    for word_length ∈ 2:6
        word = v^word_length
        Wᵢ = WhiteheadAutomorphisms(rewritingSystem, word)
        length_W = 0
        for _ in Wᵢ length_Wᵢ += 1 end
        @test length_Wᵢ == factorial(word_length) * big(2)^word_length
    end
end

function whitehead_reduce!(rewritingSystem::RewritingSystem, w::Word{T}) where {T}
    for σ ∈ WhiteheadAutomorphisms(rewritingSystem, w)
        w′ = rewrite!(one(w), w, σ)
        Base.length(w′) < Base.length(w) && return w′, σ, true
    end
    return w, nothing, false
end

function minimize!(w::Word, rewritingSystem::RewritingSystem)
    has_been_shortend_once = false
    while true
        w, _, has_been_shortened = whitehead_reduce!(rewritingSystem, w)
        has_been_shortened ? has_been_shortend_once = true : break
    end
    return w, σ, has_been_shortend_once
end


"""
Our first implementation is based on the description by Lyndon and Schupp
in their book Combinatorial Group Theory, specifically Proposition 4.19.
This version attempts to applied Whitehead transformations until the input
words no longer can be shortened.
"""
function whitehead_naive!(rewritingSystem::RewritingSystem{O, Word{T}},
                          v::Word{T},
                          w::Word{T}) where {O, T}
    # For now we assume v and w to be cyclicly reduced

    v, σv, v_has_been_shortened = minimize!(v, rewritingSystem)
    w, σw, w_has_been_shortened = minimize!(w, rewritingSystem)

    # If v and w have different lengths, there cannot exit an automorphism
    # carrying one to the other.
    Base.length(v) == Base.length(w) || return nothing
    G = automorphism_graph(A, Base.length(v))
    return path(G, v, w)
end