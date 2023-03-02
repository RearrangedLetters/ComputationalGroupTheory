#=
Models an automorphism of a free group on the given basis by prescribing images.

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

struct NielsenAutomorphisms{T}
    X::Alphabet{T}

    function NielsenAutomorphisms(X::Alphabet{T}) where {T}
        @assert issymmetric(X)
        @assert length(X) > 0
        new{T}(X)
    end
end

function Base.length(N::NielsenAutomorphisms)
    n = convert(Int, length(N.X) / 2)
    return 5n * (n - 1)
end

@enum NielsenType begin
    INVERT=1
    LEFT_MULTIPLY=2
    LEFT_MULTIPLY_INVERSE=3
    RIGHT_MULTIPLY=4
    RIGHT_MULTIPLY_INVERSE=5
end

function Base.iterate(N::NielsenAutomorphisms)
    if convert(Int, length(N.X) / 2) == 1
        return FreeGroupAutomorphism(N.X, [Word(inv(N.X, N.X[1]))]), nothing
    end
    nielsen_state = iterate(instances(NielsenType))
    return iterate(N, (1, 1, nielsen_state))
end

#=
If i == j, then nielsentype is INVERT, otherwise we have x = X[i] ≠ X[j] = y
and iterate over the NielsenTypes.
=#
function Base.iterate(N::NielsenAutomorphisms{T}, state) where {T}
    X = N.X
    n = convert(Int, length(X) / 2)
    if isnothing(state) return nothing end
    i, j, nielsenstate = state
    if i == j == n return nothing end
    if isnothing(nielsenstate)
        if j < n
            return iterate(N, (i, j + 1, iterate(instances(NielsenType))))
        elseif i < n
            return iterate(N, (i + 1, 1, iterate(instances(NielsenType))))
        else
            return nothing
        end
    end
    nielsentype, next_nielsenstate = nielsenstate
    w₁ = Vector{Word{T}}()
    for k ∈ 1:(i - 1) push!(w₁, Word(X[k])) end  # todo: shouldn't subindexing work?
    x = Word{T}(X[i])
    if i == j && nielsentype ≠ INVERT
        if j < n
            return iterate(N, (i, j + 1, nielsenstate))
        elseif i < n
            return iterate(N, (i + 1, 1, nielsenstate))
        else
            return nothing
        end
    end
    if nielsentype == INVERT
        w₂ = inv(x, X)
    else
        y = Word{T}(X[j])
        w₂ =    if nielsentype == LEFT_MULTIPLY          y * x
            elseif nielsentype == LEFT_MULTIPLY_INVERSE  inv(y, X) * x
            elseif nielsentype == RIGHT_MULTIPLY         x * y
            elseif nielsentype == RIGHT_MULTIPLY_INVERSE x * inv(y, X)
        end
    end
    w₃ = Vector{Word{T}}()
    for k ∈ (i + 1):n push!(w₃, Word(X[k])) end
    return FreeGroupAutomorphism(X, [w₁; [w₂]; w₃]),
           (i, j, iterate(instances(NielsenType), next_nielsenstate))
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
            image = if t[index] == 1 x  # todo: use @enum
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