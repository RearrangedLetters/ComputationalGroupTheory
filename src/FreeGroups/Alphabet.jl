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

function Alphabet(letters::Vector{T}, inverses::Vector{T}) where {T}
    A = Alphabet([letters; inverses])
    @assert length(letters) == length(inverses)
    for i ∈ 1:length(letters)
        setinverse!(A, letters[i], inverses[i])
    end
    return A
end

Base.getindex(A::Alphabet{T}, letter::T) where {T} = A.positions[letter]  # return ordinal of letter, i.e. A[a] -> 1
Base.getindex(A::Alphabet, index::Integer) = A.letters[index]  # return n-th letter

function Base.push!(A::Alphabet{T}, letter::T) where {T}
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
	print(io, "Alphabet = ", Tuple(A.letters))
    #= println("Alphabet with ", length(A), " letters of type ", T, ":")
    for a ∈ A
        println(a)
    end =#
end

function enumeratewords(A::Alphabet, wordlength)
    collect(Iterators.product(ntuple(_ -> A.letters, wordlength)...))
end

issymmetric(A::Alphabet) = all([hasinverse(A, l) for l ∈ A])

#=
Expects a string with all lowercase characters. Each character then is taken
as a Symbol and it's inverse is set to be the uppercase version.
This string macro is compatible with the word macro in Word.jl
=#
macro symmetric_alphabet_str(string::String)
    @assert all([islowercase(s) for s ∈ string])
    letters  = [Symbol(s) for s ∈ string]
    inverses = [Symbol(s) for s ∈ uppercase(string)]
    return :(Alphabet($letters, $inverses))
end