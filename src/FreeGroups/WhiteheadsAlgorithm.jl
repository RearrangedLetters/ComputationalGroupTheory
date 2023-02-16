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

struct WhiteheadAutomorphisms{T}
    X::Alphabet{T}
    X⁺⁻::Vector{T}
    w::Word{T}
    iterator

    function WhiteheadAutomorphisms(X::Alphabet{T}, w::Word{T}) where {T}
        X⁻ = [inv(x) for x ∈ W.X.letters]
        rank = length(X)
        # This is an iterator over all words of length rank - 1 over letters 1 through 4:
        iterator = Iterators.product(ntuple(_ -> 1:4, rank - 1)...)
        new{T}(X, [X.letters; X⁻], w, iterator)
    end
end

function iterate(W::WhiteheadAutomorphisms)
    return FreeGroupAutomorphism(W.X, W.X.letters),
           (1, iterate(powerset(W.X.letters)))
end

"""
First we iterate over the Nielsen automorphisms. This immediatedly has the consequence
that our implementation of Whitehead's algorithm employs the Nielsen-first heuristic.
Since Nielsen automorphisms are Whitehead automorphisms, these will be considered twice.
In the grand scheme of things is only quadratic additional work that won't be done in
more than 99% of cases.

The state is a tuple (i, j, n, iₙ, m, iₗ, iterator_state) consisting of:
    • i and j correspond to the i-th and j-th letter in X
    • n = 1,...,5n(n-1) counting how many Nielsen automorphisms we already considered.
      This decides, when we start considering Whitehead autormorphisms of type i.
    • iₙ ∈ {1, 2, 3, 4, 5} corresponding to one of the possible Nielsen automorphisms
      after we fixed two basis elements.
    • m - 1 is the number of
    • iₗ defines the position of the fixed element in the last loop
    • iterator_state is the state of W.iterator
"""
function iterate(W::WhiteheadAutomorphisms{T}, state) where {T}
    i, j, n, iₙ, m, iₗ, iterator_state = state
    rank = length(X)
    X = W.X
    if n ≤ 5 * rank * (rank - 1)
        return iterate(N::NielsenAutomorphisms, state)
    elseif m ≤ factorial(2 * rank)
        σ_images = nthperm(W.X⁺⁻, m)[begin:rank]
        σ = FreeGroupAutomorphism(X, σ_images)
        return σ, (i, j, n, iₙ, m + 1)
    elseif l ≤ 2 * rank * 4^(rank - 1) - 2 * rank
        σ_images = Vector{Word{T}}
        index = 1
        iteration = iterate(W.iterator, iterator_state)
        if isnothing(iteration)
            return iterate(W, (i, j, n, iₙ, m, iₗ + 1, iterate(iterator)))
        else 
            t, s = iteration
        end
        a = Word[X[iₗ]]
        while index ≤ 2 * rank
            index == iₗ && continue
            x = Word(X[index])
            image = if t[index] == 1 x
            elseif t[index] == 2 x * a
            elseif t[index] == 3 inv(a) * x
            elseif t[index] == 4 inv(a) * x * a end
            push!(σ_images, image)
        end
        σ = FreeGroupAutomorphism(X, σ_images)
        return σ, (i, j, n, iₙ, m, iₗ, s)
    end
    return nothing
end

@testset "Count Free Group Automorphisms" begin
    """
    First we assert that the iteration protocol actually does the desired number
    of iterations. This number is not equal to the actual number of Whitehead
    automorphisms.
    """
    W₁ = WhiteheadAutomorphisms(X, word"a")
    length_W₁ = 0
    for _ in W₁ length_W₁ += 1 end
    @test length_W₁ == 2

    v = word"a"

    for wordlength ∈ 2:5
        word = v^wordlength
        Wᵢ = WhiteheadAutomorphisms(rewritingSystem, word)
        length_Wᵢ = 0
        for _ in Wᵢ length_Wᵢ += 1 end
        @test length_Wᵢ == factorial(wordlength) * big(2)^wordlength
    end

    """
    Now we check if we obtain the correct number of Whitehead automorphisms.
    According to [HAR] there are 
    """
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

"""
For now see the description of iterate(::WhiteheadAutomorphisms, state)
"""
function iterate(N::NielsenAutomorphisms, state)
    i, j, n, iₙ, _ = state
    rank = length(X)
    i > rank || return
    j > rank || return
    if n ≤ 5 * rank * (rank - 1)
        j > rank && return iterate(N, (i + 1, 1, n, iₙ))
        # The case i > rank should never occur since we also count and check if n is in bounds.
        w₁ = Vector{Word{T}}()
        for k ∈ 1:(i - 1)
            push!(w₁, Word(X[k]))
        end
        x = Word{T}(X[i])
        y = Word{T}(X[j])
        w₂ =    if iₙ == 1 inv(x)
            elseif iₙ == 2 y * x
            elseif iₙ == 3 inv(y) * x
            elseif iₙ == 4 x * y
            elseif iₙ == 5 x * inv(y)
            else return iterate(N, (i, j, n, 1)) end
        w₃ = Vector{Word{T}}()
        for k ∈ (i + 1):length(N.X)
            push!(w₃, Word(N.X[k]))
        end
        σ = FreeGroupAutomorphism(X, append!(w₁, append!([w₂], w₃)))
        return σ, (i, j + 1, n + 1, mod1(iₙ + 1, 5))
    else
        return nothing
    end
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
By the definition of the iteration protocol, Nielsen autormorphisms are
considered first for minimization. This implementation thus employs the
Nielsen-first heuristic. This strategy has been shown to perform like a
polynomial time algorithm in experiments, making it practical for many
applications. The worst-case complexity however is still exponential.
The output is a list of automorphisms that need to be composed in
reverse order to get the desired automorphism. This composition can be
calculated with the method compose below. It shall however be noted, that
this composition might exhibit exponential image length.
"""
function whitehead!(X::Alphabet{T}, v::Word{T}, w::Word{T}) where {T}
    # For now we assume v and w to be cyclicly reduced

    v, _, _ = minimize!(X, v)
    w, _, _ = minimize!(X, w)

    # If v and w have different lengths, there cannot exit an automorphism
    # carrying one to the other.
    Base.length(v) == Base.length(w) || return nothing
    G = automorphism_graph(X, Base.length(v))
    return connect_depthfirst(X, v, w)
end

function whitehead(X::Alphabet{T}, v::Word{T}, w::Word{T}) where {T}
    return whitehead!(X, copy(v), copy(w))
end

function isprimitive(X::Alphabet, w::Word)
    τ = whitehead(X, w, X[1])
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

"""
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
"""