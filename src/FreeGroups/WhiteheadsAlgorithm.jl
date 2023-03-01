using ComputationalGroupTheory
using Test
using Combinatorics

#=
The type Word can be used like a cyclic word by using getcyclicindex instead
of getindex. To avoid bugs that are hard to find, the usually indexing is not
cyclic by default.
=#
begin
    w = word"abba"
    cyclicword = Vector{Symbol}()
    for i ∈ 1:8
        push!(cyclicword, getcyclicindex(w, i))
    end
    @test cyclicword == word"abbaabba"
end

function whitehead_reduce!(X::Alphabet{T}, w::Word{T}) where {T}
    for σ ∈ WhiteheadAutomorphisms(X)
        w′ = rewrite(σ(w), X)
        length(w′) < length(w) && return w′, σ, true
    end
    return w, nothing, false
end

function minimize!(X::Alphabet{T}, w::Word{T}) where {T}
    has_been_shortend_once = false
    S = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_shortened = whitehead_reduce!(X, w)
        if has_been_shortened
            has_been_shortend_once = true
            push!(S, σ)
        end
    end
    return w, S, has_been_shortend_once
end

struct AutomorphismGraph{T}
    X::Alphabet{T}
    vertices::Vector{Word{T}}
    vertex_indices::Dict{Word{T}, Int}
    edges::Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}

    function AutomorphismGraph(X::Alphabet{T};
                               wordlength::Int,
                               automorphisms=WhiteheadAutomorphisms(X)) where {T}
        vertices = Vector{Word{T}}()
        vertex_indices = Dict{Word{T}, Int}()
        edges = Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}()

        i = 1
        for w ∈ Words(X, wordlength)
            push!(vertices, w)
            push!(vertex_indices, w => i)
            i += 1
        end

        sizehint!(edges, length(vertices))
        for _ ∈ 1:length(vertices)
            push!(edges, Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}())
        end
        for v ∈ vertices
            for σ ∈ automorphisms
                t = σ(v)
                push!(edges[vertex_indices[v]], σ => t)
            end
        end

        new{T}(X, vertices, vertex_indices, edges)
    end
end

getindex(G::AutomorphismGraph{T}, w::Word{T}) where {T} = G.vertex_indices[w]
order(G::AutomorphismGraph) = length(G.vertices)
Base.size(G::AutomorphismGraph)  = sum(length.(G.edges))
edges(G::AutomorphismGraph) = G.edges
edges(G::AutomorphismGraph, s::Word{T}) where {T} = G.edges[G[s]]
wordlength(G::AutomorphismGraph) = length(first(G.vertices))
typeof(::AutomorphismGraph{T}) where {T} = T

#=
An edge (σ, w) is in the output, if σ(v) = w.
=#
function edges(G::AutomorphismGraph{T}, v::Word{T}, w::Word{T}) where {T}
    if !haskey(G.vertex_indices, v) || !haskey(G.vertex_indices, v)
        return Pair{FreeGroupAutomorphism{T}, Word{T}}[]
    end
    return filter(e -> e[2] == w, G.edges[G.vertex_indices[v]])
end

function Base.in(v::Word, G::AutomorphismGraph)
    typeof(v) ≠ typeof(G) && return false
    length(v) ≠ wordlength(G) && return false
    return v ∈ G.vertices
end

#=
Use a depth first search to find a path from s to t. The reverse composition of the
automorphisms in τ then satisfies s ↦ t. This algorithm only needs additional memory
linear in the length of the path, but potentially finds longer paths than a shortest
path algorithm. However, a shortest path algorithm needs exponential additional memory.
Whichever version has the more desirable behavior needs to be determined experimentally.

This algorithm technically modifies its inputs, however these are only auxillary data
structures and should never be passed by a user. The first three inputs aren't being
modified.
=#
function connect_depthfirst(G::AutomorphismGraph{T}, s::Word{T}, t::Word{T},
                            visited=falses(order(G)), τ=FreeGroupAutomorphism[]) where {T}
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

#=
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
=#
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

#=
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
=#