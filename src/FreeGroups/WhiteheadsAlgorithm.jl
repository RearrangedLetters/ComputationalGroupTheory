using Combinatorics

@doc """
    whitehead_reduce(w::Vector{Word}, X::Basis [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the combined length of the given words in w with the given automorphisms with
respect to the given free group basis X. If the reduction has been successful, the shortened word,
the shortening automorphism and true is returned (indicating that the new word is indeed shorter).
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a minimizing
automorphism will be found if and only if there is one.
Always returns a cyclically reduced word.
"""
function whitehead_reduce(w::Vector{Word{T}}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    startlength = sum(length.(w))
    for σ ∈ automorphisms
        v = σ.(w)
        cyclically_reduce!(v, alphabet(X))
        if sum(length.(v)) < startlength
            return v, σ, true
        end
    end
    return cyclically_reduce!(deepcopy(w), alphabet(X)), nothing, false
end

@doc """
    whitehead_reduce(w::Word, X::Basis, automorphisms [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the length of the given word w once with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found if and only if there is one.
"""
function whitehead_reduce(w::Word{T}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    v, σ, has_been_reduced = whitehead_reduce([w], X; automorphisms=automorphisms)
    return v..., σ, has_been_reduced
end

@doc """
    minimize!(w::Vector{Word{T}}, X::Basis{T} [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to minimize the combined length of the given words in w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism(s) Σ and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
Σ is a list of automorphisms whose composition will reduce w to an automorphically equivalent word
list of minimal length (under autormophisms).
"""
function minimize!(w::Vector{Word{T}}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    has_been_reduced_once = false
    Σ = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_reduced = whitehead_reduce(w, X, automorphisms=automorphisms)
        if has_been_reduced
            has_been_reduced_once = true
            push!(Σ, σ)
        else
            break
        end
    end
    return w, reverse(Σ), has_been_reduced_once
end

@doc """
    minimize!(w::Word{T}, X::Basis{T} [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to minimize the length of the given word w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism(s) Σ and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
Σ is a list of automorphisms whose composition will reduce w to an automorphically equivalent word
of minimal length (under autormophisms).
"""
function minimize!(w::Word{T}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    w, Σ, has_been_minimized = minimize!([w], X; automorphisms=automorphisms)
    return w..., Σ, has_been_minimized
end

abstract type AbstractAutomorphismGraph{T} end

vertexindex(G::AbstractAutomorphismGraph{T}, w::Vector{Word{T}}) where {T} = G.vertex_indices[w]
order(G::AbstractAutomorphismGraph) = length(G.vertices)
vertices(G::AbstractAutomorphismGraph) = G.vertices
Base.size(G::AbstractAutomorphismGraph)  = sum(length.(G.edges))
edges(G::AbstractAutomorphismGraph) = G.edges
wordlength(G::AbstractAutomorphismGraph) = length.(first(G.vertices))
typeof(::AbstractAutomorphismGraph{T}) where {T} = T

@doc """
    edges(G, v)

Return the outgoing edges of v ∈ G. The caller is responsible for calling this function with
a vertex in the graph, otherwise an error is thrown and an empty result of the correct type
is returned.
"""
function edges(G::AbstractAutomorphismGraph{T}, v) where {T}
    if haskey(G.vertex_indices, v)
        return G.edges[G.vertex_indices[v]]
    end
    @error "Vertex not in graph!"
    return Tuple{FreeGroupAutomorphism{T}, Vector{Word{T}}}[]
end

@doc """
    edges(G, v, w)

Return all edges (σ, w) ∈ G leading from v to w. An edge is in the graph iff σ(v) = w.
Possibly return multiple edges, as the automorphism graph is a multi-graph.
"""
function edges(G::AbstractAutomorphismGraph{T}, v::Vector{Word{T}}, w::Vector{Word{T}}) where {T}
    if haskey(G.vertex_indices, v)
        return filter(e -> arecyclicallyequal(e[2], w), edges(G, v))
    end
    @error "Vertex $v not in graph!"
end

@doc (@doc edges(::AbstractAutomorphismGraph{T}, ::Vector{Word{T}}, ::Vector{Word{T}}) where {T})
function edges(G::AbstractAutomorphismGraph{T}, v::Word{T}, w::Word{T}) where {T}
    return edges(G, [v], [w])
end

@doc """
    ∈(w::Word, G::AbstractAutomorphismGraph)

Return if w (understood as a cyclic word) is among the vertices of G.
"""
function Base.in(w::Vector{Word{T}}, G::AbstractAutomorphismGraph{T}) where {T}
    length.(w) ≠ wordlength(G) && return false
    for v ∈ vertices(G)
        if arecyclicallyequal(w, v) return true end
    end
    return false
end

@doc (@doc in(::Vector{Word{T}}, ::AbstractAutomorphismGraph{T}) where {T})
Base.in(w::Word{T}, G::AbstractAutomorphismGraph{T}) where {T} = [w] ∈ G

@doc """
    ∈(w::Word, G::AbstractAutomorphismGraph)

Decide if w (understood as a cyclic word) is in G.
"""
function Base.in(w::Word{T}, words::Vector{Word{T}}) where {T}
    for v ∈ words
        if arecyclicallyequal(w, v) return true end
    end
    return false
end

@doc """
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

@doc """
A multi graph with vertives representing tuples of (possibly) cyclic words and directed
edges labeled by automorphisms taking the words of one vertex componentwise to the other
vertex.
"""
struct AutomorphismGraph{T} <: AbstractAutomorphismGraph{T}
    X::Basis{T}
    vertices::Vector{Vector{Word{T}}}
    wordlengths::Vector{Int}
    vertex_indices::Dict{Vector{Word{T}}, Int}
    edges::Vector{Vector{Pair{FreeGroupAutomorphism{T}, Vector{Word{T}}}}}

    function AutomorphismGraph(X::Basis{T};
        wordlengths::Vector{Int},
        automorphisms=WhiteheadAutomorphisms(X),
        usecyclicwords=true) where {T}

        vertices = Vector{Vector{Word{T}}}()
        vertex_indices = Dict{Vector{Word{T}}, Int}()

        i = 1
        for w ∈ Words(alphabet(X), sum(wordlengths))
            v = splitbefore(w, wordlengths)
            for vᵢ ∈ v cyclically_reduce!(vᵢ, alphabet(X)) end
            if !usecyclicwords || v ∉ vertices
                if length.(v) == wordlengths
                    push!(vertices, v)
                    push!(vertex_indices, v => i)
                    i += 1
                end
            end
        end

        edges = Vector{Vector{Pair{FreeGroupAutomorphism{T}, Vector{Word{T}}}}}()
        sizehint!(edges, length(vertices))
        for _ ∈ 1:length(vertices)
            push!(edges, Vector{Pair{FreeGroupAutomorphism{T}, Word{T}}}())
        end
        for s ∈ vertices
            for σ ∈ automorphisms
                t = Vector{Word{T}}()
                for sᵢ ∈ s push!(t, cyclically_reduce!(σ(sᵢ), alphabet(X))) end
                if length.(t) == wordlengths && s ≠ t
                    push!(edges[vertex_indices[s]], σ => t)
                end
            end
        end

        new{T}(X, vertices, wordlengths, vertex_indices, edges)
    end
end

@doc """
    connect_depthfirst(G::AbstractAutomorphismGraph, s::Word, t::Word)

Use a depth first search to find a path from s to t. The reverse composition of the
automorphisms in τ then satisfies s ↦ t. This algorithm only needs additional memory
linear in the length of the path, but potentially finds longer paths than a shortest
path algorithm. However, a shortest path algorithm needs exponential additional memory.
Whichever version has the more desirable behavior needs to be determined experimentally.
"""
function connect_depthfirst(G::AbstractAutomorphismGraph{T}, s::Vector{Word{T}}, t::Vector{Word{T}}) where {T}
    if s ∉ G || t ∉ G
        return nothing
    end
    visited = falses(order(G))
    visited[vertexindex(G, s)] = true
    return connect_depthfirst!(G, s, t, visited, FreeGroupAutomorphism{T}[])
end

@doc (@doc connect_depthfirst)
function connect_depthfirst!(G::AbstractAutomorphismGraph{T}, s::Vector{Word{T}}, t::Vector{Word{T}},
                            visited::BitVector, τ::Vector{FreeGroupAutomorphism{T}}) where {T}
    s == t && return τ
    for (σ, v) ∈ edges(G, s)
        iᵥ = vertexindex(G, v)
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

@doc (@doc connect_depthfirst)
function connect_depthfirst(G::AbstractAutomorphismGraph{T}, s::Word{T}, t::Word{T}) where {T}
    return connect_depthfirst(G, [s], [t])
end

@doc (@doc whitehead_naive(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_naive!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, σ₁, _ = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))
    w, σ₂, _ = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))

    length.(v) ≠ length.(w) && return false, σ₁, σ₂, FreeGroupAutomorphism{T}[]
    G = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=WhiteheadAutomorphisms(X))
    τ = connect_depthfirst(G, v, w)
    return !isempty(τ), σ₁, σ₂, τ
end

@doc (@doc whitehead_naive(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_naive!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!([v], [w], X)
end

@doc (@doc whitehead_nielsenfirst(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenfirst!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, σ₁₁, _  = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))
    v, σ₁₂, _  = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))
    w, σ₂₁, _  = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))
    if length(w) ≠ length(v)
        w, σ₂₂, _  = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))
    else
        σ₂₂ = FreeGroupAutomorphism{T}[]
    end
    
    σ₁ = [σ₁₂; σ₁₁]
    σ₂ = [σ₂₂; σ₂₁]

    length.(v) ≠ length.(w) && return false, σ₁, σ₂, FreeGroupAutomorphism{T}[]

    G₁ = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=NielsenAutomorphisms(X))
    τ₁ = connect_depthfirst(G₁, v, w)
    if !isnothing(τ₁) return true, σ₁, σ₂, τ₁ end
    G₂ = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=WhiteheadAutomorphisms(X))
    τ₂ = connect_depthfirst(G₂, v, w)
    return !isempty(τ₂), σ₁, σ₂, τ₂
end

@doc (@doc whitehead_nielsenfirst(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenfirst!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!([v], [w], X)
end

@doc (@doc whitehead_nielsenonly(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenonly!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, σ₁, _ = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))
    w, σ₂, _ = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))

    length.(v) ≠ length.(w) && return FreeGroupAutomorphism{T}[]

    G = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=NielsenAutomorphisms(X))
    τ = connect_depthfirst(G, v, w)
    return !isempty(τ), σ₁, σ₂, τ
end

@doc (@doc whitehead_nielsenonly(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenonly!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!([v], [w], X)
end

@doc """
    whitehead_naive(v::Word, w::Word, X::Basis)

Implements Whitehead's algorithm by using all Whitehead automorphisms.

Decide, if there is an automorphism carrying v to w. If such an automorphism exists,
return (true, σ₁, σ₂, τ) such that τ(σ₁(v)) = σ₂(w). The automorphisms σ₁ & σ₂ minimize
the lengths of v and w respectively; τ maps one minimized word to the other.
If no automorphism carrying v to w exists, then (false, σ₁, σ₂, nothing) is returned.
Here again are v and w minimized by σ₁ & σ₂ respectively.

This version of Whitehead's algorithm has exponential time & memory requirements, consider
using whitehead_nielsenfirst(v, w, X) instead.
"""
function whitehead_naive(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_naive!(deepcopy(v), deepcopy(w), X)
end

@doc (@doc whitehead_naive(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T})
function whitehead_naive(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!([deepcopy(v)], [deepcopy(w)], X)
end

@doc """
    whitehead_nielsenfirst(v, w, X::Basis)

Decide, if there is an automorphism carrying v to w. If such an automorphism exists,
return (true, σ₁, σ₂, τ) such that τ(σ₁(v)) = σ₂(w). The automorphisms σ₁ & σ₂ minimize
the lengths of v and w respectively; τ maps one minimized word to the other.
If no automorphism carrying v to w exists, then (false, σ₁, σ₂, nothing) is returned.
Here again are v and w minimized by σ₁ & σ₂ respectively.

This version of Whitehead's algorithm uses the Nielsen-first heuristic. In most cases
Nielsen automorphisms are sufficient, so these are tried first. In the worst-case however,
it is still necessary to use Whitehead automorphisms. Hence this version also has
exponential time and memory requirements
"""
function whitehead_nielsenfirst(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!(deepcopy(v), deepcopy(w), X)
end

@doc (@doc whitehead_nielsenfirst(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenfirst(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!([deepcopy(v)], [deepcopy(w)], X)
end

@doc """
    whitehead_nielsenonly(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses only Nielsen automorphisms.
The result might not always be correct.
"""
function whitehead_nielsenonly(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!(deepcopy(v), deepcopy(w), X)
end

@doc (@doc whitehead_nielsenonly(::Vector{Word{T}}, ::Vector{Word{T}}, ::Basis{T}) where {T})
function whitehead_nielsenonly(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!([deepcopy(v)], [deepcopy(w)], X)
end

@doc """
    isprimitive_naive(w::Word, X)

Decide, if w is primitive in the free group with basis X using Whitehead automorphisms.
"""
function isprimitive_naive(w::Word{T}, X::Basis{T}) where {T}
    if length(w) == 1 return true elseif length(w) == 0 return false end
    v, _, _ = minimize!(deepcopy(w), X; automorphisms=WhiteheadAutomorphisms(X))
    return length(v) == 1
end

@doc """
    isprimitive_nielsenfirst(w::Word, X)

Decide, if w is primitive in the free group with basis X using Nielsen automorphisms
and then Whitehead automorphisms.
"""
function isprimitive_nielsenfirst(w::Word{T}, X::Basis{T}) where {T}
    if length(w) == 1 return true elseif length(w) == 0 return false end
    v, _, _ = minimize!(deepcopy(w), X; automorphisms=NielsenAutomorphisms(X))
    if length(v) == 1 return true end
    minimize!(v, X; automorphisms=WhiteheadAutomorphisms(X))
    return length(v) == 1
end

@doc """
    isprimitive_nielsenonly(w::Word, X)

Decide, if w is primitive in the free group with basis X using only Nielsen automorphisms.
The result might not be correct, but at least is with high probability.
"""
function isprimitive_nielsenonly(w::Word{T}, X::Basis{T}) where {T}
    if length(w) == 1 return true elseif length(w) == 0 return false end
    v, _, _ = minimize!(deepcopy(w), X; automorphisms=NielsenAutomorphisms(X))
    return length(v) == 1
end

@doc """
    isprimitive(w::Word, X)

Decide, if w is primitive in the free group with basis X using the most practical
version of Whitehead's algorithm.
"""
isprimitive = isprimitive_nielsenfirst

@doc """
    whitehead(w::Word, v::Word, X::Basis)
"""
whitehead = whitehead_nielsenfirst

"""
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
    [LYN] "Combinatorial Group Theory" by R. C. Lyndon and P. E. Schupp
    [RIC] Lecture: www.youtube.com/watch?v=dQw4w9WgXcQ
"""