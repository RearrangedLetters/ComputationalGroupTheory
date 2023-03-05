using ComputationalGroupTheory
using Test
using Combinatorics


function whitehead_reduce!(w::Word{T}, X::Basis{T}) where {T}
    for σ ∈ WhiteheadAutomorphisms(X)
        w′ = cyclically_reduce(σ(w), X.alphabet)
        length(w′) < length(w) && return w′, σ, true
    end
    return w, nothing, false
end

function minimize!(w::Word{T}, X::Basis{T}) where {T}
    has_been_shortened = false
    S = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_shortened = whitehead_reduce!(w, X)
        if has_been_shortened
            push!(S, σ)
        else
            break
        end
    end
    return w, S, has_been_shortened
end

struct AutomorphismGraph{T}
    X::Basis{T}
    vertices::Vector{Word{T}}
    vertex_indices::Dict{Word{T}, Int}
    edges::Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}

    function AutomorphismGraph(X::Basis{T};
                               wordlength::Int,
                               automorphisms=WhiteheadAutomorphisms(X)) where {T}
        vertices = Vector{Word{T}}()
        vertex_indices = Dict{Word{T}, Int}()
        
        i = 1
        for w ∈ Words(X.alphabet, wordlength)
            push!(vertices, w)
            push!(vertex_indices, w => i)
            i += 1
        end
        
        edges = Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}()
        sizehint!(edges, length(vertices))
        for _ ∈ 1:length(vertices)
            push!(edges, Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}())
        end
        for v ∈ vertices
            for σ ∈ automorphisms
                t = cyclically_reduce(σ(v), X.alphabet)
                if length(t) == wordlength
                    push!(edges[vertex_indices[v]], σ => t)
                end
            end
        end

        new{T}(X, vertices, vertex_indices, edges)
    end
end

getindex(G::AutomorphismGraph{T}, w::Word{T}) where {T} = G.vertex_indices[w]
order(G::AutomorphismGraph) = length(G.vertices)
vertices(G::AutomorphismGraph) = G.vertices
Base.size(G::AutomorphismGraph)  = sum(length.(G.edges))
edges(G::AutomorphismGraph) = G.edges
wordlength(G::AutomorphismGraph) = length(first(G.vertices))
typeof(::AutomorphismGraph{T}) where {T} = T

function edges(G::AutomorphismGraph, v::Word{T}) where {T}
    if haskey(G.vertex_indices, v)
        return G.edges[G.vertex_indices[v]]
    end
    return Pair{FreeGroupAutomorphism{T}, Word{T}}[]  # todo: this should throw an error
end

#=
An edge (σ, w) is in the output, if σ(v) = w.
=#
function edges(G::AutomorphismGraph{T}, v::Word{T}, w::Word{T}) where {T}
    if haskey(G.vertex_indices, v) || !haskey(G.vertex_indices, v)
        return filter(e -> e[2] == w, G.edges[G.vertex_indices[v]])
    end
    return Pair{FreeGroupAutomorphism{T}, Word{T}}[]   # todo: this should throw an error
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
=#
function connect_depthfirst(G::AutomorphismGraph{T}, s::Word{T}, t::Word{T}) where {T}
    visited = falses(order(G))
    visited[G.vertex_indices[s]] = true
    return connect_depthfirst!(G, s, t, visited, FreeGroupAutomorphism{T}[])
end

function connect_depthfirst!(G::AutomorphismGraph{T}, s::Word{T}, t::Word{T},
                            visited::BitVector, τ::Vector{FreeGroupAutomorphism{T}}) where {T}
    s == t && return τ
    for (σ, v) ∈ edges(G, s)
        iᵥ = G.vertex_indices[v]
        if !visited[iᵥ]
            visited[iᵥ] = true
            push!(τ, σ)
            return connect_depthfirst!(G, v, t, visited, τ)
        else
            !isempty(τ) && pop!(τ)
        end
    end
    return reverse!(τ)
end

function compose(τ::Vector{FreeGroupAutomorphism{T}}) where {T}
    if isempty(τ) 
        return FreeGroupAutomorphism{T}()
    else
        return foldr(∘, τ)
    end
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
calculated with the method compose. It shall however be noted, that
this composition might exhibit exponential image length.
=#
function whitehead_naive!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    cyclically_reduce!(v, X.alphabet)
    cyclically_reduce!(w, X.alphabet)

    v, S₁, _ = minimize!(v, X)
    w, _, _ = minimize!(w, X)

    length(v) ≠ length(w) && return nothing
    G = AutomorphismGraph(X; wordlength=length(v), automorphisms=WhiteheadAutomorphisms(X))
    return [connect_depthfirst(G, v, w); S₁]
end

function whitehead_nielsenfirst!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    cyclically_reduce!(v, X.alphabet)
    cyclically_reduce!(w, X.alphabet)

    v, S₁, _ = minimize!(v, X)
    w, _, _ = minimize!(w, X)

    length(v) ≠ length(w) && return nothing

    G₁ = AutomorphismGraph(X; wordlength=length(v), automorphisms=NielsenAutomorphisms(X))
    τ₁ = connect_depthfirst(G₁, v, w)
    if !isnothing(τ₁) return τ₁ end
    G₂ = AutomorphismGraph(X; wordlength=length(v), automorphisms=WhiteheadAutomorphisms(X))
    return [connect_depthfirst(G₂, v, w); S₁]
end

function whitehead_naive(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!(deepcopy(v), deepcopy(w), X)
end

function whitehead_nielsenfirst(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!(deepcopy(v), deepcopy(w), X)
end

function isprimitive_naive(w::Word, X::Basis)
    v = Word(X[1])
    if w == v return true end
    τ = whitehead_naive(w, v, X)
    return isnothing(τ) ? false : length(τ) > 0
end

function isprimitive_nielsenfirst(w::Word, X::Basis)
    v = Word(X[1])
    if w == v return true end
    τ = whitehead_nielsenfirst(w, v, X)
    return isnothing(τ) ? false : length(τ) > 0
end

#=
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
=#