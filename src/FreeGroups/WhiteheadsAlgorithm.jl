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

struct FreeGroupAutomorphism{T}
    basis::Alphabet{T}
    images::Vector{Word{T}}

    function FreeGroupAutomorphism(basis::Alphabet{T}, images::Vector{Word{T}}) where {T}
        @assert length(basis) == length(images)
        new{T}(basis, images)
    end

    function FreeGroupAutomorphism{T}() where {T}
        new{T}(Alphabet{T}(), Vector{Word{T}}())
    end
end

function apply!(σ::FreeGroupAutomorphism{T}, w::Word{T}) where {T}
    return replace_all!(w, σ.basis, σ.images)
end

function push!(σ::FreeGroupAutomorphism{T}, replacement::Pair{T, Word{T}})
    letter, word = replacement
    push!(σ.basis, letter)
    push!(σ.images, word)
    return σ
end

(σ::FreeGroupAutomorphism{T})(w::Word{T}) where {T} = apply!(σ, deepcopy(w))


struct WhiteheadAutomorphisms
    rewritingSystem::RewritingSystem
    w::Word
    wordlength::Int
    number_of_permutations::BigInt
    number_of_subsets::BigInt

    function WhiteheadAutomorphisms(rewritingSystem::RewritingSystem, w::Word)
        wordlength = Base.length(w)
        new(rewritingSystem, w, wordlength, factorial(big(wordlength)), big(2)^wordlength)
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

    for wordlength ∈ 2:6
        word = v^wordlength
        Wᵢ = WhiteheadAutomorphisms(rewritingSystem, word)
        length_Wᵢ = 0
        for _ in Wᵢ length_Wᵢ += 1 end
        @test length_Wᵢ == factorial(wordlength) * big(2)^wordlength
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

struct NielsenAutomorphisms
    rewritingSystem::RewritingSystem
    w::Word

    function NielsenAutomorphisms(rewritingSystem::RewritingSystem, w::Word)
        new(rewritingSystem, w)
    end
end

function iterate(N::NielsenAutomorphisms, state)
    i, j = state
    if i ≤ length(N.w)
        x = N.w[i]
        v = Word([x])
        for y ∈ N.w
            y ≠ x || break
            j ≤ 5 || return iterate(N, (i + 1, 1))

            w = Word(if j == 1 [inv(x)]
                 elseif j == 2 [y, x]
                 elseif j == 3 [inv(y), x]
                 elseif j == 4 [x, y]
                 elseif j == 5 [x, inv(y)]
                 end)
            return w
        end
    end
    return nothing
end

@testset "Count Nielsen Automorphisms" begin
    v = word"a"

    for wordlength ∈ 1:6
        word = v^wordlength
        Nᵢ = NielsenAutomorphisms(rewritingSystem, word)
        length_Nᵢ = 0
        for _ in Nᵢ length_Nᵢ += 1 end
        @test length_Nᵢ == 5 * wordlength^2
    end
end

struct AutomorphismGraph{T}
    A::Alphabet{T}
    vertices::Vector{Word{T}}
    vertex_indices::Dict{Word{T}, Int}
    edges::Vector{Pair(FreeGroupAutomorphism{T}, Vector{Word{T}})}

    function AutomorphismGraph{T}(A::Alphabet{T}, wordlength::Int) where {T}
        numvertices = big(length(A))^wordlength
        resize!(vertices, numvertices)
        resize!(vertex_indices, numvertices)
        resize!(edges, numvertices)

        i = 1
        for w ∈ enumeratewords(A, wordlength)
            push!(vertices, Word(collect(w)))
            push!(vertex_indices, (w, i))
            i += 1
        end

        for v ∈ vertices
            for σ ∈ WhiteheadAutomorphisms(rws, v)
                t = σ(v)
                iₜ = findfirst(x -> x == t, vertices)
                push!(edges[iₜ], σ => t)
            end
        end
    end
end

getindex(G::AutomorphismGraph{T}, w::Word{T}) where {T} = G.vertex_indices[w]
order(G::AutomorphismGraph) = length(G.vertices)
edges(G::AutomorphismGraph) = G.edges
edges(G::AutomorphismGraph, s::Word{T}) = G.edges[G[s]]

"""
Use a depth first search to find a path from s to t. The reverse composition of the
automorphisms in τ then satisfies s ↦ t. This algorithm only needs additional memory linear in the
length of the path, but potentially finds longer paths than a shortest path algorithm.
However, a shortest path algorithm needs exponential additional memory. Whichever
version has the more desirable behavior needs to be determined experimentally.

This algorithm technically modifies its inputs, however these are only auxillary data
structures and should never be passed by a user. The first three inputs aren't being
modified.
"""
function connect_depthfirst(G::AutomorphismGraph{T}, s::Word{T}, t::Word{T},
                            visited=falses(order(G)), τ=FreeGroupAutomorphism[])
    s == t && return τ

    for (σ, v) ∈ edges(G, s)
        iᵥ = G[w]
        if !visited[iᵥ]
            visited[iᵥ] = true
            push!(τ, σ)
            return connect_depthfirst(G, v, t, visited, τ)
        end
    end

    return τ
end

"""
Our first implementation is based on the description by Lyndon and Schupp
in their book Combinatorial Group Theory, specifically Proposition 4.19.
This version attempts to applied Whitehead transformations until the input
words no longer can be shortened. Then a path in the automorphism graph, if
it exists, defines a desired automorphism.
The output is actually a list of automorphisms that need to be applied in
reverse order.
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
    return connect_depthfirst(G, v, w)
end

function whitehead_naive(rewritingSystem::RewritingSystem{O, Word{T}},
                         v::Word{T},
                         w::Word{T}) where {O, T}
    return whitehead_naive!(rewritingSystem, copy(v), copy(w))
end

function isirreducible_naive(X::Alphabet, w::Word)
    return length(whitehead_naive(X, w, X[1])) > 0
end