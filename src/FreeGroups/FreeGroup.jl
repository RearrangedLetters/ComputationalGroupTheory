abstract type AbstractFinitelyPresentedGroup <: Group end
abstract type GroupElement end

struct Basis{T} <: AbstractAlphabet{T}
    alphabet::Alphabet{T}

    function Basis(A::Alphabet{T}) where {T}
        @assert issymmetric(A)
        n = convert(Int, length(A) / 2)
        @assert all(i -> inv(A.letters[i], A) == A.letters[n + i], 1:n)
        new{T}(A)
    end
end

Base.length(basis::Basis) = convert(Int, length(basis.alphabet) / 2)
Base.getindex(basis::Basis, i) = basis.alphabet[i]
alphabet(basis::Basis) = basis.alphabet
letters(basis::Basis) = basis.alphabet.letters

abstract type AbstractRelation{T} end

struct Relation{T} <: AbstractRelation{T}
	relation::Pair{Word{T}, Word{T}}
end

struct FinitelyPresentedGroup{T} <: AbstractFinitelyPresentedGroup
	basis::Basis{T}
	relations::Vector{<:AbstractRelation{T}}
end