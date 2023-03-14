abstract type AbstractAlphabet{T} end

struct Alphabet{T} <: AbstractAlphabet{T}
	letters::Vector{T}
	positions::Dict{T, Int}
    inverses::Dict{T, T}  # or store a vector with placeholders and keep these consistent

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
    @assert length(inverses) ≤ length(letters)
    A = Alphabet([letters; inverses])
    for i ∈ 1:length(inverses)
        setinverse!(A, letters[i], inverses[i])
    end
    return A
end

Base.getindex(A::Alphabet{T}, letter::T) where {T} = A.positions[letter]  # return ordinal of letter, i.e. A[a] -> 1
Base.getindex(A::Alphabet, index::Integer) = A.letters[index]  # return n-th letter

"""
    push!(A::Alphabet{T}, letter::T)

Add a letter to the alphabet.
"""
function Base.push!(A::Alphabet{T}, letter::T) where {T}
    push!(A.letters, letter)
    position[A.letters[last]] = lastindex(A.letters)
    return A
end

"""
    setinverse!(A::Alphabet{T}, x::T, X::T)

Set the inverse of x ∈ A to be X ∈ A and vice versa. If this is called
with one letter not in A, then the function crashes.
"""
function setinverse!(A::Alphabet{T}, x::T, X::T) where {T}
    A.inverses[x] = X
    A.inverses[X] = x
end

Base.inv(letter::T, A::Alphabet{T}) where {T} = A.inverses[letter]
Base.inv(index::Integer, A::Alphabet{T}) where {T} = A.inverses[letters[index]]

hasinverse(A::Alphabet{T}, letter::T) where {T} = haskey(A.inverses, letter)
hasinverse(A::Alphabet, index::Integer) = hasinverse(A, A[index])

function isinverse(A::Alphabet{T}, letter₁::T, letter₂::T) where {T}
    if hasinverse(A, letter₁) && hasinverse(A, letter₂)
        return inv(letter₁, A) == letter₂
    else
        return false
    end
end

function Base.iterate(A::Alphabet, i=1)
    return if i > length(A) nothing else (A.letters[i], i + 1) end
end

Base.length(A::Alphabet) = length(A.letters)

function Base.show(io::IO, A::Alphabet{T}, oneliner=false) where {T}
    if oneliner
	    print(io, "Alphabet = ", Tuple(A.letters))
    else
        println("Alphabet with ", length(A), " letters of type ", T, ":")
        for a ∈ A
            println(a)
        end
    end
end

issymmetric(A::Alphabet) = all([hasinverse(A, l) for l ∈ A])

"""
    @alphabet"..."

Turn each individual character in the string into a Symbol and create an alphabet with
those letters. Duplicate letters will appear only once.
"""
macro alphabet_str(string::String)
    return :(Alphabet([Symbol(s) for s ∈ Set($string)]))
end

"""
    @symmetric_alphabet(string)

Expects a string with unique, all lowercase characters. Each character then is taken
as a Symbol and it's inverse is set to be the uppercase form of the character.
This string macro is compatible with the word macro in Word.jl

# Examples
@symmetric_alphabet"ab" creates an instance of Alphabet{Symbol} with letters
:a, :b, :A, :B in this exact order.
"""
macro symmetric_alphabet_str(string::String)
    @assert all([islowercase(s) for s ∈ string])
    letters  = [Symbol(s) for s ∈ string]
    inverses = [Symbol(s) for s ∈ uppercase(string)]
    return :(Alphabet($letters, $inverses))
end