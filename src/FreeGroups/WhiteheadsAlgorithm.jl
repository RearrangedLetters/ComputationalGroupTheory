include("Word.jl")
include("Alphabet.jl")
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

#=
Models an automorphism of a free group on the given basis.

Struct-invariant:
    • length(images) == length(basis) / 2
=#
struct FreeGroupAutomorphism{T}
    basis::Alphabet{T}
    images::Vector{Word{T}}

    function FreeGroupAutomorphism(basis::Alphabet{T}, images::Vector{Word{T}}) where {T}
        if length(images) == convert(Int, length(basis) / 2)
            images = [images; [inv(w, basis) for w ∈ images]]
        end
        new{T}(basis, images)
    end

    function FreeGroupAutomorphism(basis::Alphabet{T}) where {T}
        new{T}(basis, [Word(x) for x ∈ basis])
    end

    function FreeGroupAutomorphism{T}() where {T}
        new{T}(Alphabet{T}(), Vector{Word{T}}())
    end
end

basis(σ::FreeGroupAutomorphism) = σ.basis
images(σ::FreeGroupAutomorphism) = σ.images

function Base.:(==)(σ::FreeGroupAutomorphism, τ::FreeGroupAutomorphism)
    return basis(σ) == basis(τ) && images(σ) == images(τ)
end

function apply!(σ::FreeGroupAutomorphism{T}, w::Word{T}) where {T}
    return replace_all!(w, basis(σ).letters[1:length(images(σ))], images(σ))
end

(σ::FreeGroupAutomorphism{T})(w::Word{T}) where {T} = apply!(σ, Base.copy(w))
(σ::FreeGroupAutomorphism{T})(x::T) where {T} = apply!(σ, Word(x))

function Base.push!(σ::FreeGroupAutomorphism{T}, replacement::Pair{T, Word{T}}) where {T}
    letter, word = replacement
    push!(basis(σ), letter)
    push!(images(σ), word)
    return σ
end

function Base.show(io::IO, σ::FreeGroupAutomorphism)
    print("Free group automorphism on ", σ.basis, " with mapping: ")
    for i ∈ 1:length(σ.images)
        print(io, σ.basis[i], " ↦ ", σ.images[i], ", ")
    end
end

struct WhiteheadAutomorphisms{T}
    X::Alphabet{T}

    function WhiteheadAutomorphisms(X::Alphabet{T}) where {T}
        @assert issymmetric(X)
        new{T}(X)
    end
end

#=
Constructs the non-trivial Nielsen automorphisms where x = X[i] is mapped to
one of the five (indicated by iₙ) possible values where y = X[j].
=#
function nielsen(A::Alphabet{T}, i::Int, j::Int, iₙ::Int) where {T}
    if i == j @assert iₙ == 1 end
    # if iₙ == 1 @assert i == j end
    X = A.letters[1:convert(Int, length(A) / 2)]
    w₁ = Vector{Word{T}}()
    for k ∈ 1:(i - 1) push!(w₁, Word(X[k])) end  # todo: use copy
    x = Word{T}(X[i])
    if iₙ == 1
        w₂ = inv(x, A)
    else
        y = Word{T}(X[j])
        w₂ =    if iₙ == 2 y * x
            elseif iₙ == 3 inv(y, A) * x
            elseif iₙ == 4 x * y
            elseif iₙ == 5 x * inv(y, A)
            else return nothing end
    end
    w₃ = Vector{Word{T}}()
    for k ∈ (i + 1):length(X) push!(w₃, Word(X[k])) end  # todo: use copy
    σ_images = [w₁; [w₂]; w₃]
    σ = FreeGroupAutomorphism(A, σ_images)
    return σ
end

function Base.iterate(W::WhiteheadAutomorphisms{T}) where {T}
    # This is an iterator over all words of length (rank - 1) over letters 1 through 4:
    iterator = Iterators.product(ntuple(_ -> 1:4, length(W.X) - 1)...)
    _, iterator_state = Base.iterate(iterator)
    return FreeGroupAutomorphism(W.X),  # the intial call gives the identity
           (i=1, j=1, n=1, iₙ=1, m=1, l=1, iₗ=1, iterator, iterator_state)
end

#=
First we iterate over the Nielsen automorphisms. This immediatedly has the consequence
that our implementation of Whitehead's algorithm employs the Nielsen-first heuristic.
Since Nielsen automorphisms are Whitehead automorphisms, these will be considered twice.
In the grand scheme of things this amounts to only quadratic additional work that won't
be done in more than 99% of cases [HAR].

The state is a tuple (i, j, n, iₙ, m, iₗ, iterator_state) consisting of:
    • i and j correspond to the i-th and j-th letter in X
    • n = 1,...,5n(n-1) counting how many Nielsen automorphisms we already considered.
      This decides, when we start considering Whitehead autormorphisms of type i.
    • iₙ ∈ {1, 2, 3, 4, 5} corresponding to one of the possible Nielsen automorphisms
      after we fixed two basis elements.
    • m - 1 counts the number of Whitehead automorphisms of the first type considered so far
    • l - 1 counts the number of Whitehead automorphisms of the second type considered so far
    • iₗ defines the position of the fixed element in the last loop
    • iterator is an Iterator obtained from Base.iterate(::WhiteheadAutomorphisms) (see above)
    • iterator_state is the state of iterator
=#
function Base.iterate(W::WhiteheadAutomorphisms{T}, state) where {T}
    i, j, n, iₙ, m, l, iₗ, iterator, iterator_state = state
    X = W.X
    rank = convert(Int, length(X) / 2)

    # First we cover the Nielsen automorphisms
    if n ≤ 5 * rank * (rank - 1) || n ≤ 1
        j > rank && return iterate(W, (i + 1, 1, n, 1, m, l, iₗ, iterator, iterator_state))
        i == j && iₙ ≠ 1 && return iterate(W, (i, j + 1, n, iₙ, m, l, iₗ, iterator, iterator_state))
        σ = nielsen(X, i, j, iₙ)
        if isnothing(σ)
            return iterate(W, (i, j + 1, n, 1, m, l, iₗ, iterator, iterator_state))
        else 
            return σ, (i, j, n + 1, iₙ + 1, 1, m, l, iₗ, iterator, iterator_state)
        end
    
    # Then we cover the Whitehead autormorphisms of the permutation type
    elseif m ≤ factorial(2 * rank)
        m == 1 && return iterate(W, (i, j, n, iₙ, m + 1, l, iₗ, iterator, iterator_state))  # Skip identity
        σ_image_vector = nthperm(X.letters, m)
        σ_images = [Word(letter) for letter ∈ σ_image_vector]
        σ = FreeGroupAutomorphism(X, σ_images)
        return σ, (i, j, n, iₙ, m + 1, l, iₗ, iterator, iterator_state)

    # Lastly we cover the Whitehead automorphisms of the multiplication type
    elseif l ≤ 2 * rank * 4^(rank - 1) - 2 * rank
        σ_images = Vector{Word{T}}()
        index = 1
        iteration = iterate(iterator, iterator_state)
        if isnothing(iteration)
            return iterate(W, (i, j, n, iₙ, m, l, iₗ + 1, iterator, Base.iterate(iterator)))
        else 
            t, _ = iteration
        end
        a = Word(X[iₗ])
        while index ≤ 2 * rank
            index == iₗ && continue
            x = Word(X[index])
            image = if t[index] == 1 x
                elseif t[index] == 2 x * a
                elseif t[index] == 3 inv(a, X) * x
                elseif t[index] == 4 inv(a, X) * x * a end
            push!(σ_images, image)
        end
        σ = FreeGroupAutomorphism(X, σ_images)
        return σ, (i, j, n, iₙ, m, l + 1, iₗ, iterate, iterator_state)
    end
    return nothing
end

function Base.length(W::WhiteheadAutomorphisms)
    length_W::BigInt = 0
    for _ in W
        length_W += 1
    end
    return length_W
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
    edges::Vector{Pair{FreeGroupAutomorphism{T}, Vector{Word{T}}}}

    function AutomorphismGraph(X::Alphabet{T}, wordlength::Int) where {T}
        numvertices = big(length(X))^wordlength
        # Nonsense with push!:
        #= resize!(vertices, numvertices)
        resize!(vertex_indices, numvertices)
        resize!(edges, numvertices) =#

        i = 1
        for w ∈ enumeratewords(X, wordlength)
            push!(vertices, Word(collect(w)))
            push!(vertex_indices, (w, i))
            i += 1
        end

        for v ∈ vertices
            for σ ∈ WhiteheadAutomorphisms(X)
                t = σ(v)
                iₜ = findfirst(x -> x == t, vertices)
                push!(edges[iₜ], σ => t)
            end
        end

        new{T}(X, vertices)
    end
end

getindex(G::AutomorphismGraph{T}, w::Word{T}) where {T} = G.vertex_indices[w]
order(G::AutomorphismGraph) = length(G.vertices)
edges(G::AutomorphismGraph) = G.edges
edges(G::AutomorphismGraph, s::Word{T}) where {T} = G.edges[G[s]]

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