using ComputationalGroupTheory
using Test
using Combinatorics

"""
    whitehead_reduce(w::Word{T}, X::Basis{T}, automorphisms [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the length of the given word w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
"""
function whitehead_reduce(w::Word{T}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    w′ = copy(w)
    for σ ∈ automorphisms
        w′ = cyclically_reduce(σ(w′), alphabet(X))
        length(w′) < length(w) && return w′, σ, true
    end
    return w, nothing, false
end

"""
    tuple_whitehead_reduce(w::Vector{Word{T}}, X::Basis{T} [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the combined length of the given words in w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a minimizing
automorphism will be found if and only iff there is one.
"""
function tuple_whitehead_reduce(w::Vector{Word{T}}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    w′ = copy(w)
    for σ ∈ automorphisms
        w′ = (v -> cyclically_reduce(v, alphabet(X))).(w′)
        length(sum(length.(w′))) < sum(length.(w)) && return w′, σ, true  # todo: is this the correct check?
    end
    return w, nothing, false
end

"""
    minimize!(w::Word{T}, X::Basis{T} [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to minimize the length of the given word w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism(s) S and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
S is a list of automorphisms whose composition will reduce w to an automorphically equivalent word
of minimal length (under autormophisms).
"""
function minimize!(w::Word{T}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    has_been_shortened = false
    S = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_shortened = whitehead_reduce(w, X, automorphisms=automorphisms)
        if has_been_shortened
            push!(S, σ)
        else
            return w, S, has_been_shortened
        end
    end
    return w, nothing, false
end

"""
    minimize!(w::Vector{Word{T}}, X::Basis{T} [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to minimize the combined length of the given words in w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism(s) S and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
S is a list of automorphisms whose composition will reduce w to an automorphically equivalent word
list of minimal length (under autormophisms).
"""
function tuple_minimize!(w::Vector{Word{T}}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    has_been_shortened = false
    S = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_shortened = tuple_whitehead_reduce(w, X, automorphisms=automorphisms)
        if has_been_shortened
            push!(S, σ)
        else
            return w, S, has_been_shortened
        end
    end
    return w, nothing, false
end

"""
    AbstractAutomorphismGraph{T}
Define an interface to a graph 
"""
abstract type AbstractAutomorphismGraph{T} end

getindex(G::AbstractAutomorphismGraph{T}, w::Word{T}) where {T} = G.vertex_indices[w]
order(G::AbstractAutomorphismGraph) = length(G.vertices)
vertices(G::AbstractAutomorphismGraph) = G.vertices
Base.size(G::AbstractAutomorphismGraph)  = sum(length.(G.edges))
edges(G::AbstractAutomorphismGraph) = G.edges
wordlength(G::AbstractAutomorphismGraph) = length(first(G.vertices))
typeof(::AbstractAutomorphismGraph{T}) where {T} = T

"""
    edges(G, v)

Return the outgoing edges from v ∈ G.

"""
function edges(G::AbstractAutomorphismGraph, v::Word{T}) where {T}
    if haskey(G.vertex_indices, v)
        return G.edges[G.vertex_indices[v]]
    end
    @error "Vertex not in graph!"
end

"""
    edges(G, v, w)

Return all edges (σ, w) ∈ G leading from v to w. An edge is in the graph, iff σ(v) = w.
Possibly return multiple edges as the automorphism graph is a multi-graph.
"""
function edges(G::AbstractAutomorphismGraph{T}, v::Word{T}, w::Word{T}) where {T}
    if haskey(G.vertex_indices, v)
        return filter(e -> arecyclicallyequal(e[2], w), edges(G, v))
    end
    @error "Vertex $v not in graph!"
end

"""
    ∈(w::Word, G::AbstractAutomorphismGraph)

Return if w (understood as a cyclic word) is among the vertices of G.
"""
function Base.in(w::Word{T}, G::AbstractAutomorphismGraph{T}) where {T}
    length(w) ≠ wordlength(G) && return false

    for v ∈ G.vertices
        if arecyclicallyequal(w, v) return true end
    end
    return false
end

"""
    ∈(w::Word, G::AbstractAutomorphismGraph)

Decide if w (understood as a cyclic word) is in G.
"""
function Base.in(w::Word{T}, words::Vector{Word{T}}) where {T}
    for v ∈ words
        if arecyclicallyequal(w, v) return true end
    end
    return false
end

"""
    ∈(edge, edges)

Check if the given edge is among the given edges.
"""
function Base.in(edge::Pair{FreeGroupAutomorphism{T}, Word{T}}, 
                 edges::Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}) where {T}
    (σ, w) = edge
    for (τ, v) ∈ edges
        if σ == τ && arecyclicallyequal(w, v) return true end
    end
    return false
end

"""
    SimpleAutomorphismGraph

Represent a graph where the vertices are words in the free group with
basis X and each edge is labeled by an automorphism that takes the origin
word to the image under this automorphism.
It is possible to take only the cylic words as vertices.
"""
struct SimpleAutomorphismGraph{T} <: AbstractAutomorphismGraph{T}
    X::Basis{T}
    vertices::Vector{Word{T}}
    vertex_indices::Dict{Word{T}, Int}
    edges::Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}

    """
        SimpleAutomorphismGraph(X::Basis, wordlength, automorphisms[; usecyclicwords=true])

    Construct a graph with vertices either being all words or all cyclic words of the
    given word length; the directed edges are labeled by the automorphisms taking one
    word to another.
    """
    function SimpleAutomorphismGraph(X::Basis{T};
                               wordlength::Int,
                               automorphisms=WhiteheadAutomorphisms(X);
                               usecyclicwords=true) where {T}
        vertices = Vector{Word{T}}()
        vertex_indices = Dict{Word{T}, Int}()
        
        i = 1
        for w ∈ Words(alphabet(X), wordlength)
            if !usecyclicwords || w ∉ vertices
                push!(vertices, w)
                push!(vertex_indices, w => i)
                i += 1
            end
        end
        
        edges = Vector{Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}}()
        sizehint!(edges, length(vertices))
        for _ ∈ 1:length(vertices)
            push!(edges, Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}())
        end
        for v ∈ vertices
            for σ ∈ automorphisms
                t = cyclically_reduce(σ(v), alphabet(X))
                if length(t) == wordlength
                    push!(edges[vertex_indices[v]], σ => t)
                end
            end
        end

        new{T}(X, vertices, vertex_indices, edges)
    end
end

"""

"""
struct AutomorphismGraph{T} <: AbstractAutomorphismGraph{T}
    X::Basis{T}
    vertices::Vector{Vector{Word{T}}}
    wordlengths::Vector{Int}
    vertex_indices::Dict{Vector{Word{T}}, Int}
    edges::Vector{Vector{Pair{FreeGroupAutomorphism{T}, Vector{Word{T}}}}}

    function AutomorphismGraph(X::Basis{T};
        wordlengths::Vector{Int},
        automorphisms=WhiteheadAutomorphisms(X)) where {T}

        vertices = Vector{Vector{Word{T}}}
        vertex_indices = Dict{Word{T}, Int}()

        i = 1
        for w ∈ Words(alphabet(X), sum(wordlength))
            v = splitbefore(w, wordlengths)
            for vᵢ ∈ v cyclically_reduce!(vᵢ, alphabet(X)) end
            if !usecyclicwords || w ∉ vertices
                push!(vertices, w)
                push!(vertex_indices, w => i)
                i += 1
            end
        end

        edges = Vector{Vector{Pair{FreeGroupAutomorphism{T}, Vector{Word{T}}}}}()
        sizehint!(edges, length(vertices))
        for _ ∈ 1:length(vertices)
            push!(edges, Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}())
        end
        for v ∈ vertices
            for σ ∈ automorphisms
                for vᵢ ∈ v cyclically_reduce!(σ(vᵢ), alphabet(X)) end
                if length.(v) == wordlengths
                    push!(edges[vertex_indices[v]], σ => t)
                end
            end
        end

        new{T}(X, vertices, wordlengths, vertex_indices, edges)
    end
end

"""
    connect_depthfirst(G::AbstractAutomorphismGraph, s::Word, t::Word)

Use a depth first search to find a path from s to t. The reverse composition of the
automorphisms in τ then satisfies s ↦ t. This algorithm only needs additional memory
linear in the length of the path, but potentially finds longer paths than a shortest
path algorithm. However, a shortest path algorithm needs exponential additional memory.
Whichever version has the more desirable behavior needs to be determined experimentally.
"""
function connect_depthfirst(G::AbstractAutomorphismGraph{T}, s::Word{T}, t::Word{T}) where {T}
    visited = falses(order(G))
    visited[G.vertex_indices[s]] = true
    return connect_depthfirst!(G, s, t, visited, FreeGroupAutomorphism{T}[])
end

function connect_depthfirst!(G::AbstractAutomorphismGraph{T}, s::Word{T}, t::Word{T},
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

"""
    compose(τ::Vector{FreeGroupAutomorphism})

Compose the list of automorphisms returned by connect_depthfirst into a single
automorphism.
"""
function compose(τ::Vector{FreeGroupAutomorphism{T}}) where {T}
    if isempty(τ) 
        return FreeGroupAutomorphism{T}()
    else
        return foldr(∘, τ)
    end
end

function whitehead_naive!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, S₁, _ = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))

    length(v) ≠ length(w) && return nothing
    G = SimpleAutomorphismGraph(X; wordlength=length(v), automorphisms=WhiteheadAutomorphisms(X))
    return [connect_depthfirst(G, v, w); S₁]
end

"""

"""
function whitehead_nielsenfirst!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, _, _  = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))
    v, S₁, _ = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))

    length(v) ≠ length(w) && return nothing

    G₁ = SimpleAutomorphismGraph(X; wordlength=length(v), automorphisms=NielsenAutomorphisms(X))
    τ₁ = connect_depthfirst(G₁, v, w)
    if !isnothing(τ₁) return τ₁ end
    G₂ = SimpleAutomorphismGraph(X; wordlength=length(v), automorphisms=WhiteheadAutomorphisms(X))
    return [connect_depthfirst(G₂, v, w); S₁]
end

"""

"""
function whitehead_nielsenonly!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, S₁, _ = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))

    length(v) ≠ length(w) && return nothing

    G₁ = SimpleAutomorphismGraph(X; wordlength=length(v), automorphisms=NielsenAutomorphisms(X))
    τ₁ = connect_depthfirst(G₁, v, w)
    return [τ₁; S₁]
end

"""
    whitehead_naive(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses the Whitehead autormophisms.
"""
function whitehead_naive(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!(deepcopy(v), deepcopy(w), X)
end

"""
    whitehead_nielsenfirst(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses Nielsen automorphisms first
and only then uses the Whitehead autormophisms.
"""
function whitehead_nielsenfirst(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!(deepcopy(v), deepcopy(w), X)
end

"""
    whitehead_nielsenonly(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses only Nielsen automorphisms.
The result might not always be correct. todo: Check this and maybe @warn
"""
function whitehead_nielsenonly(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!(deepcopy(v), deepcopy(w), X)
end

"""
    isprimitive_naive(w::Word, X)

Decide, if w is primitive in the free group with basis X using Whitehead automorphisms.
"""
function isprimitive_naive(w::Word{T}, X::Basis{T}) where {T}
    v = Word(X[1])
    if w == v return true end
    τ = whitehead_naive(w, v, X)
    return isnothing(τ) ? false : length(τ) > 0
end

"""
    isprimitive_nielsenfirst(w::Word, X)

Decide, if w is primitive in the free group with basis X using Nielsen automorphisms
and then Whitehead automorphisms.
"""
function isprimitive_nielsenfirst(w::Word{T}, X::Basis{T}) where {T}
    v = Word(X[1])
    if w == v return true end
    τ = whitehead_nielsenfirst(w, v, X)
    return isnothing(τ) ? false : length(τ) > 0
end

"""
    isprimitive_nielsenonly(w::Word, X)

Decide, if w is primitive in the free group with basis X using only Nielsen automorphisms.
The result might not be correct, but at least is with high probability.
"""
function isprimitive_nielsenonly(w::Word{T}, X::Basis{T}) where {T}
    @warn "The result of isprimitive_nielsenonly is not always correct."
    v = Word(X[1])
    if w == v return true end
    τ = whitehead_nielsenfirst(w, v, X)
    return isnothing(τ) ? false : length(τ) > 0
end

"""
    isprimitive(w::Word, X)

Decide, if w is primitive in the free group with basis X using the most practical
version of Whitehead's algorithm
"""
isprimitive = isprimitive_nielsenfirst

"""
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
    [LYN] "Combinatorial Group Theory" by R. C. Lyndon and P. E. Schupp
    [RIC] Lecture: www.youtube.com/watch?v=dQw4w9WgXcQ
"""