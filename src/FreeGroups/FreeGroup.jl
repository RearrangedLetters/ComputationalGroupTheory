abstract type AbstractFinitelyPresentedGroup <: Group end

struct Basis{T} <: AbstractAlphabet{T}
    alphabet::Alphabet{T}

    function Basis(A::Alphabet{T}) where {T}
        @assert issymmetric(A)
        n = convert(Int, length(A) / 2)
        @assert all(i -> inv(A.letters[i], A) == A.letters[n + i], 1:n)
        new{T}(A)
    end
end

alphabet(basis::Basis) = basis.alphabet
Base.getindex(basis::Basis, i) = alphabet(basis)[i]
Base.length(basis::Basis) = convert(Int, length(alphabet(basis)) / 2)
letters(basis::Basis) = basis.alphabet.letters
Base.inv(x::T, basis::Basis{T}) where {T} = inv(x, alphabet(basis))
Base.inv(w::Word{T}, basis::Basis{T}) where {T} = inv(w, alphabet(basis))

abstract type AbstractRelation{T} end

struct Relation{T} <: AbstractRelation{T}
	relation::Pair{Word{T}, Word{T}}
end

struct FinitelyPresentedGroup{T} <: AbstractFinitelyPresentedGroup
	basis::Basis{T}
	relations::Vector{<:AbstractRelation{T}}
end