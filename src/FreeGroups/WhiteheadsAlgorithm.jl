using Combinatorics

"""
    whitehead_reduce(w::Vector{Word}, X::Basis [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the combined length of the given words in w with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a minimizing
automorphism will be found if and only iff there is one.
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
    return w, nothing, false
end

"""
    whitehead_reduce(w::Word, X::Basis, automorphisms [; automorphisms=WhiteheadAutomorphisms(X)])

Attempt to reduce the length of the given word w once with the automorphisms with respect
to the given free group basis X. If the reduction has been successful, the shortened word, the
shortening automorphism and true is returned, indicating that the new word is indeed shorter.
Otherwise, the original word, nothing and false are returned as no automorphism could be found.
By default, the Whitehead autormorphisms are being used, as they guarentee that a word of minimal
length will be found.
"""
function whitehead_reduce(w::Word{T}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    v, σ, has_been_shortened = whitehead_reduce([w], X; automorphisms=automorphisms)
    return v..., σ, has_been_shortened
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
function minimize!(w::Vector{Word{T}}, X::Basis{T}; automorphisms=WhiteheadAutomorphisms(X)) where {T}
    has_been_shortened = false
    S = Vector{FreeGroupAutomorphism{T}}()
    while true
        w, σ, has_been_shortened = whitehead_reduce(w, X, automorphisms=automorphisms)
        if has_been_shortened
            push!(S, σ)
        else
            break
        end
    end
    return w, S, has_been_shortened
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
    w, S, has_been_minimized = minimize!([w], X; automorphisms=automorphisms)
    return w..., S, has_been_minimized
end

"""
    AbstractAutomorphismGraph{T}
Define an interface to a graph 
"""
abstract type AbstractAutomorphismGraph{T} end

vertexindex(G::AbstractAutomorphismGraph{T}, w::Vector{Word{T}}) where {T} = G.vertex_indices[w]
order(G::AbstractAutomorphismGraph) = length(G.vertices)
vertices(G::AbstractAutomorphismGraph) = G.vertices
Base.size(G::AbstractAutomorphismGraph)  = sum(length.(G.edges))
edges(G::AbstractAutomorphismGraph) = G.edges
wordlength(G::AbstractAutomorphismGraph) = length.(first(G.vertices))
typeof(::AbstractAutomorphismGraph{T}) where {T} = T

"""
    edges(G, v)

Return the outgoing edges from v ∈ G.

"""
function edges(G::AbstractAutomorphismGraph, v)
    if haskey(G.vertex_indices, v)
        return G.edges[G.vertex_indices[v]]
    end
    @error "Vertex not in graph!"
end

"""
    edges(G, v, w)

Return all edges (σ, w) ∈ G leading from v to w. An edge is in the graph, iff σ(v) = w.
Possibly return multiple edges, as the automorphism graph is a multi-graph.
"""
function edges(G::AbstractAutomorphismGraph{T}, v::Vector{Word{T}}, w::Vector{Word{T}}) where {T}
    if haskey(G.vertex_indices, v)
        return filter(e -> arecyclicallyequal(e[2], w), edges(G, v))
    end
    @error "Vertex $v not in graph!"
end

function edges(G::AbstractAutomorphismGraph{T}, v::Word{T}, w::Word{T}) where {T}
    return edges(G, [v], [w])
end

"""
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

Base.in(w::Word{T}, G::AbstractAutomorphismGraph{T}) where {T} = [w] ∈ G

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

"""
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

function connect_depthfirst(G::AbstractAutomorphismGraph{T}, s::Word{T}, t::Word{T}) where {T}
    return connect_depthfirst(G, [s], [t])
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

function whitehead_naive!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, S, _ = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))

    length.(v) ≠ length.(w) && return FreeGroupAutomorphism{T}[]
    G = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=WhiteheadAutomorphisms(X))
    τ = connect_depthfirst(G, v, w)
    if isnothing(τ)
        return false, S
    else
        return true, [τ; S]
    end
end

function whitehead_naive!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!([v], [w], X)
end

"""

"""
function whitehead_nielsenfirst!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, _, _  = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))  # todo: this S is also important for the result
    v, S, _  = minimize!(v, X, automorphisms=WhiteheadAutomorphisms(X))  # todo: in all these methods we should also return a minimizing auto for the second word
    w, _, _  = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))
    if length(w) ≠ length(v)
        w, _, _  = minimize!(w, X, automorphisms=WhiteheadAutomorphisms(X))
    end

    length.(v) ≠ length.(w) && return false, FreeGroupAutomorphism{T}[]

    G₁ = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=NielsenAutomorphisms(X))
    τ₁ = connect_depthfirst(G₁, v, w)
    if !isnothing(τ₁) return true, [τ₁; S] end
    G₂ = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=WhiteheadAutomorphisms(X))
    τ₂ = connect_depthfirst(G₂, v, w)
    if isnothing(τ₂)
        return false, FreeGroupAutomorphism{T}[]
    else
        return true, [τ₂; S]
    end
end

function whitehead_nielsenfirst!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!([v], [w], X)
end

"""

"""
function whitehead_nielsenonly!(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    cyclically_reduce!(v, alphabet(X))
    cyclically_reduce!(w, alphabet(X))

    v, S, _ = minimize!(v, X, automorphisms=NielsenAutomorphisms(X))
    w, _,  _ = minimize!(w, X, automorphisms=NielsenAutomorphisms(X))

    length.(v) ≠ length.(w) && return FreeGroupAutomorphism{T}[]

    G = AutomorphismGraph(X; wordlengths=length.(v), automorphisms=NielsenAutomorphisms(X))
    τ = connect_depthfirst(G, v, w)
    if isnothing(τ)
        return false, S
    else
        return true, [τ; S]
    end
end

function whitehead_nielsenonly!(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!([v], [w], X)
end

"""
    whitehead_naive(v::Word, w::Word, X::Basis)

Attempt to find an automorphism τ of the free group with basis X, which maps
the word v to w. The minimization strategy uses the Whitehead autormophisms.
Return a list of these automorphisms such that their composition τ takes v to w,
or nothing if no such τ exists.
"""
function whitehead_naive(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_naive!(deepcopy(v), deepcopy(w), X)
end

function whitehead_naive(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_naive!([deepcopy(v)], [deepcopy(w)], X)
end

"""
    whitehead_nielsenfirst(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses Nielsen automorphisms first
and only then uses the Whitehead autormophisms.
Return a list of these automorphisms such that their composition τ takes v to w,
or nothing if no such τ exists.
"""
function whitehead_nielsenfirst(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!(deepcopy(v), deepcopy(w), X)
end

function whitehead_nielsenfirst(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenfirst!([deepcopy(v)], [deepcopy(w)], X)
end

"""
    whitehead_nielsenonly(v::Word, w::Word, X::Basis)

Attempt to find an automorphism of the free group with basis X, which maps
the word v to w. The minimization strategy uses only Nielsen automorphisms.
The result might not always be correct. todo: Check this and maybe @warn
"""
function whitehead_nielsenonly(v::Vector{Word{T}}, w::Vector{Word{T}}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!(deepcopy(v), deepcopy(w), X)
end

function whitehead_nielsenonly(v::Word{T}, w::Word{T}, X::Basis{T}) where {T}
    return whitehead_nielsenonly!([deepcopy(v)], [deepcopy(w)], X)
end

"""
    isprimitive_naive(w::Word, X)

Decide, if w is primitive in the free group with basis X using Whitehead automorphisms.
"""
function isprimitive_naive(w::Word{T}, X::Basis{T}) where {T}
    if length(w) == 1 return true elseif length(w) == 0 return false end
    v, _, _ = minimize!(deepcopy(w), X; automorphisms=WhiteheadAutomorphisms(X))
    return length(v) == 1
end

"""
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

"""
    isprimitive_nielsenonly(w::Word, X)

Decide, if w is primitive in the free group with basis X using only Nielsen automorphisms.
The result might not be correct, but at least is with high probability.
"""
function isprimitive_nielsenonly(w::Word{T}, X::Basis{T}) where {T}
    @warn "The result of isprimitive_nielsenonly is not always correct."
    if length(w) == 1 return true elseif length(w) == 0 return false end
    v, _, _ = minimize!(deepcopy(w), X; automorphisms=NielsenAutomorphisms(X))
    return length(v) == 1
end

"""
    isprimitive(w::Word, X)

Decide, if w is primitive in the free group with basis X using the most practical
version of Whitehead's algorithm.
"""
isprimitive = isprimitive_nielsenfirst

whitehead = whitehead_nielsenfirst

"""
The algorithms and insights are mostly based on the following sources:
    [HAR] "Heuristics for the Whitehead Minimization Problem" by Haralick et. al.
    [LYN] "Combinatorial Group Theory" by R. C. Lyndon and P. E. Schupp
    [RIC] Lecture: www.youtube.com/watch?v=dQw4w9WgXcQ
"""