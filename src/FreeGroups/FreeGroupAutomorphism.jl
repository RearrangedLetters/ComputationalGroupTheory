"""
    AbstractFreeGroupAutomorphism{T}

Represents an interface for free group automorphisms. Each implementation should
implement the following functions
    • basis(<:AbstractFreeGroupAutomorphism) which is the ordered list of letters on which the automorphism acts
    • images(<:AbstractFreeGroupAutomorphism), the image of the i-th element in the basis
    • (σ<:AbstractFreeGroupAutomorphism)(x::T) evaluates to σ(x)
    • (σ<:AbstractFreeGroupAutomorphism)(w::Word{T}) evaluates to σ(w)
"""
abstract type AbstractFreeGroupAutomorphism{T} end

"""
    FreeGroupAutomorphism

Represents a free group automorphism by prescribing images to the basis elements.
"""
struct FreeGroupAutomorphism{T} <: AbstractFreeGroupAutomorphism{T}
    basis::Basis{T}
    images::Vector{Word{T}}

    function FreeGroupAutomorphism(basis::Basis{T}, images::Vector{Word{T}}) where {T}
        @assert length(images) == length(basis)
        new{T}(basis, [images; [inv(w, basis) for w ∈ images]])
    end

    function FreeGroupAutomorphism(alphabet::Alphabet{T}) where {T}
        @assert issymmetric(alphabet)
        new{T}(alphabet, [[Word(x) for x ∈ alphabet]; Word(inv(x, alphabet)) for x ∈ alphabet])
    end

    function FreeGroupAutomorphism{T}() where {T}
        new{T}(Alphabet{T}(), Vector{Word{T}}())
    end
end

basis(σ::FreeGroupAutomorphism) = σ.basis
alphabet(σ::FreeGroupAutomorphism) = σ.basis.alphabet
images(σ::FreeGroupAutomorphism) = σ.images

function Base.:(==)(σ::FreeGroupAutomorphism, τ::FreeGroupAutomorphism)
    return basis(σ) == basis(τ) && images(σ) == images(τ)
end

function apply!(σ::FreeGroupAutomorphism{T}, w::Word{T}) where {T}
    toreplace = σ.basis.alphabet.letters
    replacements = [images(σ); [inv(y, σ.basis.alphabet) for y ∈ images(σ)]]
    return replace_all!(w, toreplace, replacements)
end

(σ::FreeGroupAutomorphism{T})(w::Word{T}) where {T} = apply!(σ, Base.copy(w))
(σ::FreeGroupAutomorphism{T})(x::T) where {T} = apply!(σ, Word(x))

"""
    inv(σ::FreeGroupAutomorphism)

Return σ⁻¹, the inverse automorphism.
"""
function Base.inv(σ::FreeGroupAutomorphism{T}) where {T}
    images = Vector{Word{T}}()
    for i ∈ 1:length(basis(σ))
        push!(images, inv(σ(X[i]), alphabet(σ)))
    end
    return FreeGroupAutomorphism(basis(σ), images)
end

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

struct NielsenAutomorphisms{T}
    X::Basis{T}

    function NielsenAutomorphisms(X::Basis{T}) where {T}
        @assert length(X) > 0
        new{T}(X)
    end
end

function Base.length(N::NielsenAutomorphisms)
    n = length(N.X)
    return 5n * (n - 1)
end

basis(N::NielsenAutomorphisms) = N.X
alphabet(N::NielsenAutomorphisms) = alphabet(N.X)

"""
See nielsen_automorphism.
"""
@enum NielsenType begin
    N_INVERT
    N_LEFT_MULTIPLY
    N_LEFT_MULTIPLY_INVERSE
    N_RIGHT_MULTIPLY
    N_RIGHT_MULTIPLY_INVERSE
end

"""
    nielsen_automorphism(X::Basis, type::NielsenType, x[, y=x])

Construct a FreeGroupAutomorphism that is a Nielsen automorphism of the given type.
The result fixes all y ∈ X, y ≠ x and, depending on the type, sends x to
    • x⁻¹   (N_INVERT),
    • yx    (N_LEFT_MULTIPLY),
    • y⁻¹x  (N_LEFT_MULTIPLY_INVERSE),
    • xy    (N_RIGHT_MULTIPLY),
    • xy⁻¹  (N_RIGHT_MULTIPLY_INVERSE).
"""
function nielsen_automorphism(X::Basis{T}, type::NielsenType, x::T, y::T=x) where {T}
    @assert type ≠ N_INVERT || x == y
    w₁ = Vector{Word{T}}()
    for k ∈ 1:(i - 1) push!(w₁, Word(X[k])) end
    x′ = Word(x)
    if type == N_INVERT
        w₂ = inv(x′, X)
    else
        y′ = Word(y)
        w₂ =    if nielsentype == N_LEFT_MULTIPLY          y′ * x′
            elseif nielsentype == N_LEFT_MULTIPLY_INVERSE  inv(y′, X) * x′
            elseif nielsentype == N_RIGHT_MULTIPLY         x′ * y′
            elseif nielsentype == N_RIGHT_MULTIPLY_INVERSE x′ * inv(y′, X)
        end
    end
    w₃ = Vector{Word{T}}()
    for k ∈ (i + 1):n push!(w₃, Word(X[k])) end
    return FreeGroupAutomorphism(X, [w₁; [w₂]; w₃])
end

"""
    iterate(::NielsenAutomorphisms)

Iterate over all possible Nielsen automorphisms.
"""
function Base.iterate(N::NielsenAutomorphisms)
    if length(basis(N)) == 1
        return FreeGroupAutomorphism(basis(N), [Word(inv(alphabet(N), alphabet(N)[1]))]), nothing
    end
    nielsen_state = iterate(instances(NielsenType))
    return iterate(N, (1, 1, nielsen_state))
end

"""
If i == j, then nielsentype is INVERT, otherwise we have x = X[i] ≠ X[j] = y
and iterate over the NielsenTypes.
"""
function Base.iterate(N::NielsenAutomorphisms{T}, state) where {T}
    X = basis(N)
    n = length(X)
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
    if i == j && nielsentype ≠ N_INVERT
        if j < n
            return iterate(N, (i, j + 1, nielsenstate))
        elseif i < n
            return iterate(N, (i + 1, 1, nielsenstate))
        else
            return nothing
        end
    end
    return nielsen_automorphism(X, nielsentype, X[i], X[j]),
           (i, j, iterate(instances(NielsenType), next_nielsenstate))
end

abstract type AbstractWhiteheadAutomorphisms{T} end

alphabet(W::AbstractWhiteheadAutomorphisms) = W.X.alphabet
basis(W::AbstractWhiteheadAutomorphisms) = W.X
rank(W::AbstractWhiteheadAutomorphisms) = length(basis(W))

struct WhiteheadAutomorphismsTypeI{T} <: AbstractWhiteheadAutomorphisms{T}
    X::Basis{T}
end

Base.length(W::WhiteheadAutomorphismsTypeI) = factorial(big(rank(W))) * 2^rank(W)

"""
    whiteheadI_automorphism(X::Basis, k<:Integer, inversions::Vector{Bool})
"""
function whiteheadI_automorphism(X::Basis, k<:Integer, inversions::Vector{Bool})
    @assert k ≤ factorial(big(length(X)))
    @assert length(X) == length(inversions)

    images = nthperm(letters(X)[1:rank(W)], k)
    for i ∈ 1:length(images)
        if inversions[i] images[i] = inv(alphabet(W), images[i]) end
    end
    return FreeGroupAutomorphism(X, [Word(x) for x ∈ images])
end

function Base.iterate(W::WhiteheadAutomorphismsTypeI)
    permutation = 1
    inversion_iterator = iterate(Words(Alphabet([true, false]), rank(W)))
    return iterate(W, (permutation, inversion_iterator))
end

function Base.iterate(W::WhiteheadAutomorphismsTypeI, state)
    permutation, inversion_iterator = state
    if isnothing(inversion_iterator)
        return iterate(W::WhiteheadAutomorphismsTypeI,
                      (permutation + 1, iterate(Words(Alphabet([true, false]), rank(W)))))
    else
        if permutation ≤ factorial(big(rank(W)))
            inversion, inversion_state = inversion_iterator
            return whiteheadI_automorphism(basis(W), permutation, inversion),
                   (permutation, iterate(Words(Alphabet([true, false]), rank(W)), inversion_state))
        else
            return nothing
        end
    end
end

struct WhiteheadAutomorphismsTypeII{T} <: AbstractWhiteheadAutomorphisms{T}
    X::Basis{T}
    type_iterator

    function WhiteheadAutomorphismsTypeII(X::Basis{T}) where {T}
        whitehead_types = collect(instances(WhiteheadType))
        new{T}(X, Words(Alphabet(whitehead_types), length(X) - 1))
    end
end

"""
    length(::WhiteheadAutomorphisms)

Calculate the total number of WhiteheadAutomorphisms of the second type based on
the formula 2n * 4ⁿ⁻¹ - 2n, where n is the rank.
"""
function Base.length(W::WhiteheadAutomorphismsTypeII)
    n = rank(W)
    return 2n * 4^(n-1) - 2n
end

@enum WhiteheadType begin
    W_IDENTITY
    W_RIGHT_MULTIPLY
    W_LEFT_MULTIPLY_INVERSE
    W_CONJUGATE
end

function Base.iterate(W::WhiteheadAutomorphismsTypeII)
    if rank(W) == 1
        return nothing
    else
        return iterate(W, (1, iterate(W.type_iterator)))
    end
end

function whiteheadII_automorphism(X::Basis{T}, multiplier_index, type_word) where {T}
    a = Word(X[multiplier_index])
    offset = 0
    images = Vector{Word{T}}()
    A = alphabet(X)
    for i ∈ 1:length(X)
        x = Word(X[i])
        if i == multiplier_index || i == multiplier_index - length(X)
            offset = 1
            push!(images, x)
        else
            whiteheadtype = type_word[i - offset]
            w =     if whiteheadtype == W_IDENTITY              x
                elseif whiteheadtype == W_RIGHT_MULTIPLY        x * a
                elseif whiteheadtype == W_LEFT_MULTIPLY_INVERSE inv(a, A) * x
                elseif whiteheadtype == W_CONJUGATE             inv(a, A) * x * a end
            push!(images, w)
        end
    end
    return FreeGroupAutomorphism(X, images)
end

function Base.iterate(W::WhiteheadAutomorphismsTypeII, state)
    isnothing(state) && return nothing
    multiplier_index, whiteheadtype_state = state
    multiplier_index > 2rank(W) && return nothing
    if isnothing(whiteheadtype_state)
        return iterate(W, (multiplier_index + 1, iterate(W.type_iterator)))
    else
        type_word, next_state = whiteheadtype_state
        i₁ = max(1, mod1(multiplier_index, rank(W)) - 1)
        i₂ = min(length(type_word), mod1(multiplier_index, rank(W)) + 1)
        if type_word[begin:i₁] == [W_IDENTITY for _ ∈ 1:i₁] &&
           type_word[i₂:end] == [W_IDENTITY for _ ∈ i₂:length(type_word)]
            return iterate(W, (multiplier_index, iterate(W.type_iterator, next_state)))
        else
            return whiteheadII_automorphism(W.X, multiplier_index, type_word),
                   (multiplier_index, iterate(W.type_iterator, next_state))
        end
    end
end

struct WhiteheadAutomorphisms{T} <: AbstractWhiteheadAutomorphisms{T}
    X::Basis{T}

    function WhiteheadAutomorphisms(X::Basis{T}) where {T}
        @assert issymmetric(X.alphabet)
        new{T}(X)
    end
end

function Base.length(W::WhiteheadAutomorphisms)
    return length(WhiteheadAutomorphismsTypeI(W.X)) +
           length(WhiteheadAutomorphismsTypeII(W.X))
end

function Base.iterate(W::WhiteheadAutomorphisms)
    W₁ = WhiteheadAutomorphismsTypeI(basis(W))
    W₂ = WhiteheadAutomorphismsTypeII(basis(W))
    state = (W₁, iterate(W₁), W₂, iterate(W₂))
    return iterate(W, state)
end

function Base.iterate(::WhiteheadAutomorphisms, state)
    isnothing(state) && return nothing
    W₁, state₁, W₂, state₂ = state
    if !isnothing(state₁)
        return state₁[1], (W₁, iterate(W₁, state₁[2]), W₂, state₂)
    elseif !isnothing(state₂)
        return state₂[1], (W₁, nothing, W₂, iterate(W₂, state₂[2]))
    else
        return nothing
    end
end

