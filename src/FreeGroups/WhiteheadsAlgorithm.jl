include("Word.jl")
include("Alphabet.jl")
using Test
using Combinatorics

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

basis(σ::FreeGroupAutomorphism) = σ.basis
images(σ::FreeGroupAutomorphism) = σ.images

function apply!(σ::FreeGroupAutomorphism{T}, w::Word{T}) where {T}
    return replace_all!(w, basis(σ), images(σ))
end

(σ::FreeGroupAutomorphism{T})(w::Word{T}) where {T} = apply!(σ, deepcopy(w))

function push!(σ::FreeGroupAutomorphism{T}, replacement::Pair{T, Word{T}})
    letter, word = replacement
    push!(basis(σ), letter)
    push!(images(σ), word)
    return σ
end

struct WhiteheadAutomorphisms
    X::Alphabet
    w::Word
    wordlength::Int
    number_of_permutations::BigInt

    function WhiteheadAutomorphisms(X::Alphabet, w::Word)
        wordlength = Base.length(w)
        new(X, w, wordlength, factorial(big(wordlength)))
    end
end

function iterate(W::WhiteheadAutomorphisms)
    return FreeGroupAutomorphism(X, letters(X)), (1, iterate(powerset(W.w.letters)))
end

function iterate(W::WhiteheadAutomorphisms, state)
    X = alphabet(W.rewritingSystem)
    i, powerset, j = state
    powerset_iterator = iterate(powerset, j)
    if i ≤ W.number_of_permutations
        σ_images = nthperm(letters(X), i)
        if !isnothing(powerset_iterator)
            subset_indices, _ = powerset_iterator
            for index ∈ subset_indices
                σ_images[index] = inv(σ_images[index])
            end
            return FreeGroupAutomorphism(X, σ_images)
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

function whitehead_reduce!(X::Alphabet{T}, w::Word{T}) where {T}
    for σ ∈ WhiteheadAutomorphisms(X, w)
        w′ = rewrite(σ(w), X)
        Base.length(w′) < Base.length(w) && return w′, σ, true
    end
    return w, nothing, false
end

function minimize!(X::Alphabet{T}, w::Word{T}) where {T}
    has_been_shortend_once = false
    σ = FreeGroupAutomorphism{T}()
    while true
        w, _, has_been_shortened = whitehead_reduce!(X, w)
        has_been_shortened ? has_been_shortend_once = true : break
    end
    return w, σ, has_been_shortend_once
end

struct NielsenAutomorphisms{T}
    X::Alphabet{T}
    w::Word{T}

    function NielsenAutomorphisms(X::Alphabet{T}, w::Word{T}) where {T}
        new(X, w)
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
    X::Alphabet{T}
    vertices::Vector{Word{T}}
    vertex_indices::Dict{Word{T}, Int}
    edges::Vector{Pair(FreeGroupAutomorphism{T}, Vector{Word{T}})}

    function AutomorphismGraph{T}(X::Alphabet{T}, wordlength::Int) where {T}
        numvertices = big(length(X))^wordlength
        resize!(vertices, numvertices)
        resize!(vertex_indices, numvertices)
        resize!(edges, numvertices)

        i = 1
        for w ∈ enumeratewords(X, wordlength)
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
function whitehead_naive!(X::Alphabet{T}, v::Word{T}, w::Word{T}) where {T}
    # For now we assume v and w to be cyclicly reduced

    v, _, _ = minimize!(X, v)
    w, _, _ = minimize!(X, w)

    # If v and w have different lengths, there cannot exit an automorphism
    # carrying one to the other.
    Base.length(v) == Base.length(w) || return nothing
    G = automorphism_graph(X, Base.length(v))
    return connect_depthfirst(X, v, w)
end

function whitehead_naive(X::Alphabet{T}, v::Word{T}, w::Word{T}) where {T}
    return whitehead_naive!(X, copy(v), copy(w))
end

function isprimitive_naive(X::Alphabet, w::Word)
    τ = whitehead_naive(X, w, X[1])
    return isnothing(τ) ? false : length(τ) > 0
end

@testset "Primitive elements in ℤ" begin
    X = Alphabet(:𝟙)
    setinverse!(X, :𝟙, :𝟙⁻)

    @test freeRewriteBV!(word"𝟙𝟙⁻", X) == word""

    """
    Now F₁ ≅ ⟨𝟙⟩ ≅ ℤ, and there are exactly two Irreducible elements,
    namely 𝟙 ≙ 1 and -1 ≙ 𝟙⁻.
    """
    @test isprimitive_naive(X, word"𝟙")
    @test isprimitive_naive(X, word"𝟙⁻")

    """
    Now we assert that neither 2 ≙ 𝟙𝟙 nor -2 ≙ 𝟙⁻𝟙⁻ are irreducible.
    """
    @test !isprimitive_naive(X, word"𝟙𝟙")
    @test !isprimitive_naive(X, word"𝟙⁻𝟙⁻")
end