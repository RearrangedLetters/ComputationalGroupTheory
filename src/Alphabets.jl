struct Alphabet{T}
	letters::Vector{T}  # or alternatively an ordered set
	positions::Dict{T, Int}
    inverses::Dict{T, T}

    function Alphabet(letters::Vector{T}) where T
        positions::Dict{T, Int}()
        for i ∈ 1:length(letters)
            positions[letters[i]] = i
        end
        new{T}(letters, positions, Dict{T, T}())
    end
end

# for implementing the error have a look at transversal where NotInOrbit Exception is defined

Base.getindex(A::Alphabet{T}, letter::T) where T = A.positions[letter]  # return ordinal of letter, i.e. A[a] -> 1
Base.getindex(A::Alphabet, index::Integer) = A.letters[index]  # return n-th letter

setinverse!(A::Alphabet{T}, x::T, X::T) where T = A.inverses[x] = X  # set value of "inv" involution

Base.inv(A::Alphabet{T}, letter::T) where T = A.inverses[letter]  # the inverse of letter as T
Base.inv(A::Alphabet{T}, index::Integer) where T = A.inverses[letters[index]]  # the ordinal of the inverse of 'n'-th letter

hasinverse(A::Alphabet{T}, letter::T) where T = hasinverse(A, A[letter])  # 
hasinverse(A::Alphabet, index::Integer) = haskey(A.inverses, A.letters[index])  # is the partially defined "inv" defined for this particular index?

Base.iterate(A::Alphabet, i=1) = A.letters[i]
Base.length(A) = length(A.letters)

function Base.show(io::IO, A::Alphabet{T}) where T
	println(io, "Alphabet of $T with $(length(A)) letters:")
	for letter in A
		print(io, A[letter], "\t, letter")
		if hasinverse(A, letter)
			# println(io, " with inverse", A[inv(A, A[letter]]))
		# else
			# println(io, "")
		end
	end
end