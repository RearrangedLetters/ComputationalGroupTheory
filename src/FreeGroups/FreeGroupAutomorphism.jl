
struct Basis{T}  # todo <: AbstractAlphabet
    alphabet::Alphabet{T}

    function Basis(A::Alphabet{T}) where {T}
        @assert issymmetric(A)
        n = convert(Int, length(A) / 2)
        @assert all(i -> inv(A, A.letters[i]) == A.letters[n + i], 1:n)
        new{T}(A)
    end
end

Base.length(X::Basis) = convert(Int, length(X.alphabet) / 2)
Base.getindex(X::Basis, i) = X.alphabet[i]
alphabet(X::Basis) = X.alphabet

struct FreeGroupAutomorphism{T}
    basis::Basis{T}
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

    function FreeGroupAutomorphism(basis::Basis{T}, images::Vector{Word{T}}) where {T}
        @assert length(basis) == length(images)
        new{T}(basis, images)
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

@enum NielsenType begin
    N_INVERT
    N_LEFT_MULTIPLY
    N_LEFT_MULTIPLY_INVERSE
    N_RIGHT_MULTIPLY
    N_RIGHT_MULTIPLY_INVERSE
end

function Base.iterate(N::NielsenAutomorphisms)
    if length(N.X) == 1
        return FreeGroupAutomorphism(N.X, [Word(inv(N.X.alphabet, N.X.alphabet[1]))]), nothing
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
    w₁ = Vector{Word{T}}()
    for k ∈ 1:(i - 1) push!(w₁, Word(X[k])) end  # todo: shouldn't subindexing work?
    x = Word(X[i])
    if i == j && nielsentype ≠ N_INVERT
        if j < n
            return iterate(N, (i, j + 1, nielsenstate))
        elseif i < n
            return iterate(N, (i + 1, 1, nielsenstate))
        else
            return nothing
        end
    end
    if nielsentype == N_INVERT
        w₂ = inv(x, X)
    else
        y = Word(X[j])
        w₂ =    if nielsentype == N_LEFT_MULTIPLY          y * x
            elseif nielsentype == N_LEFT_MULTIPLY_INVERSE  inv(y, X) * x
            elseif nielsentype == N_RIGHT_MULTIPLY         x * y
            elseif nielsentype == N_RIGHT_MULTIPLY_INVERSE x * inv(y, X)
        end
    end
    w₃ = Vector{Word{T}}()
    for k ∈ (i + 1):n push!(w₃, Word(X[k])) end
    return FreeGroupAutomorphism(X, [w₁; [w₂]; w₃]),
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
            images = nthperm(W.X.alphabet.letters[1:rank(W)], permutation)
            for i ∈ 1:length(images)
                if inversion[i] images[i] = inv(alphabet(W), images[i]) end
            end
            return FreeGroupAutomorphism(W.X, [Word(x) for x ∈ images]),
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

function construct_whiteheadII(X::Basis{T}, multiplier_index, type_word) where {T}
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
            return construct_whiteheadII(W.X, multiplier_index, type_word),
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

