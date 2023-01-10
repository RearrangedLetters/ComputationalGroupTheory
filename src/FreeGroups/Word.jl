abstract type AbstractWord{T} <: AbstractVector{T} end

one(w::AbstractWord) = one(typeof(w))
isone(w::AbstractWord) = iszero(length(w))

function Base.:*(w::AbstractWord, v::AbstractWord)
	return mul!(one(w), w, v)  # this is an interface function!
end

# resize!
# append!

# concrete implementation

struct Word{T} <: AbstractWord{T}  # <: AbstractVector{T}
	letter_indices::Vector{T}
end

Base.inv(w::AbstractWord{T}, A::Alphabet) = inv!(similar(w), w, A)

function inv!(out::AbstractWord, w::AbstractWord, A::Alphabet)
	resize!(out, length(w))
	# for letter in reverse(w) allocate vector containing reversed w
	for (i, letter) in enumerate(Iterators.reverse(w))
		out[i] = inv(letter, A)
	end
	return out
end

# Implement abstract Vector interface
Base.size(w::Word) = size(w.letters)
Base.getindex(w::Word, i::Integer) = w.letters[i]
Base.setindex!(w::Word, value, i::Int) = w.letters[i] = value

function mul!(out::AbstractWord, w::AbstractWord, v::AbstractWord)
	@assert out !== w  # out and w occupy different places in memory (actually out.letters and w.letters, not a problem because the structs are immutable)
	# resize!(out, length(w) + length(v))
	# this is now correct but doesn't allow us to do
	# mul!(a, a, b) override a with content of a * b
	out = resize!(out, 0)
	out = append!(out, w)  # this should work as appending of w::AbstractVector to out.letters::Vector{T} is defined
	out = append!(out, v)
	return out
end

function freeRewriteV1(w::Word)
    i = 1
	@inbounds while i < length(w)
		if iszero(w[i] * w[i + 1])
			w = w[begin:(i - 1)] * w[(i + 2):end]  # this part is sub-optimal
			i = max(1, i - 1)
		end
		i += 1
	end
end

function freeRewriteV2(w::Word)
	i = 1
	l = length(w)
	marked = trues(l)
	@inbounds while i < l
		if !marked[i]
			i += 1
			continue
		end
		if iszero(w[i] * w[i + 1])
			marked[i:(i + 1)] .= false
		end
		while !marked[i] && i ≥ 1
			i -= 1
		end
	end
	return 
end

# Compare performance on words like x¹⁰⁰X¹⁰⁰, x⁵⁰X¹⁰⁰x⁵⁰, (x¹⁰X¹⁰)¹⁰