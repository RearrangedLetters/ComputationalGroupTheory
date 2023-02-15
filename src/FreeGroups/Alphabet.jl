struct Alphabet{T}
	letters::Vector{T}  # or alternatively an ordered set
	positions::Dict{T, Int}
    inverses::Dict{T, T}  # or store a vector with placeholders and keep these consistent
    # todo: in any case, there should be consistency checks in the constructor

    function Alphabet(letters::Vector{T}) where {T}
        positions = Dict{T, Int}(letters[i] => i for i ∈ 1:length(letters))
        new{T}(letters, positions, Dict{T, T}())
    end

    function Alphabet(letters::T...) where {T}
        positions = Dict{T, Int}(letters[i] => i for i ∈ 1:length(letters))
        new{T}(collect(letters), positions, Dict{T, T}())
    end

    function Alphabet{T}() where {T}
        new{T}(Vector{T}(), Dict{T, Int}(), Dict{T, T}())
    end
end

# for implementing the error have a look at transversal where NotInOrbit Exception is defined

Base.getindex(A::Alphabet{T}, letter::T) where {T} = A.positions[letter]  # return ordinal of letter, i.e. A[a] -> 1
Base.getindex(A::Alphabet, index::Integer) = A.letters[index]  # return n-th letter

function push!(A::Alphabet{T}, letter::T) where {T}
    push!(A.letters)
    position[A.letters[last]] = lastindex(A.letters)
    return A
end

function setinverse!(A::Alphabet{T}, x::T, X::T) where {T}
    A.inverses[x] = X
    A.inverses[X] = x
end

Base.inv(A::Alphabet{T}, letter::T) where {T} = A.inverses[letter]  # the inverse of letter as T
Base.inv(A::Alphabet{T}, index::Integer) where {T} = A.inverses[letters[index]]  # the ordinal of the inverse of 'n'-th letter

hasinverse(A::Alphabet{T}, letter::T) where {T} = haskey(A.inverses, letter)  # 
hasinverse(A::Alphabet, index::Integer) = hasinverse(A, A[index])  # is the partially defined "inv" defined for this particular index?

function Base.iterate(A::Alphabet, i=1)
    return if i > length(A) nothing else (A.letters[i], i + 1) end
end

Base.length(A) = length(A.letters)

function Base.show(io::IO, A::Alphabet{T}) where {T}
	print(io, "Σ = ", Tuple(A.letters))
end

enumeratewords(A::Alphabet, wordlength)
    collect(Iterators.product(ntuple(_ -> A.letters, wordlength)...))
end