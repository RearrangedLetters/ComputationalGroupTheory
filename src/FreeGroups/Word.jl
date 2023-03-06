using ComputationalGroupTheory

abstract type AbstractWord{T} <: AbstractVector{T} end

Base.one(w::AbstractWord) = one(typeof(w))
isone(w::AbstractWord) = iszero(length(w))

function Base.:*(w::AbstractWord, v::AbstractWord)
	return mul!(one(w), w, v)  # this is an interface function
end

mutable struct Word{T} <: AbstractWord{T}  # todo: mutable is probably not necessary here
	letters::Vector{T}

	function Word(letters::Vector{T}) where {T}
		new{T}(letters)
	end

	function Word(w::Word{T}) where {T}
		new{T}(deepcopy(w.letters))
	end

	function Word{T}(w::Word{T}) where {T}
		new{T}(deepcopy(w.letters))
	end

	function Word(letter::T) where {T}
		new{T}([letter])
	end

	function Word(letters::T...) where {T}
		new{T}(collect(letters))
	end

	function Word{T}() where {T}
		new{T}(Vector{T}())
	end
end

function Base.resize!(w::Word, size::Integer)
	Base.resize!(w.letters, size)
	return w
end

function Base.similar(w::AbstractWord, ::Type, dims::Base.Dims{1})
	ans = one(w)
	resize!(ans, first(dims))
	return ans
end

function Base.append!(w::Word{T}, v::Word{T}) where {T}
	Base.append!(w.letters, v.letters)
	return w
end

Base.one(::Word{T}) where {T} = Word{T}()
Base.inv(w::AbstractWord{T}, A::Alphabet) where {T} = inv!(similar(w), w, A)

function inv!(out::AbstractWord, w::AbstractWord, A::Alphabet)
	resize!(out, length(w))
	# for letter in reverse(w) allocate vector containing reversed w
	for (i, letter) in enumerate(Iterators.reverse(w))
		out[i] = inv(A, letter)
	end
	return out
end

# Implement abstract Vector interface
Base.size(w::Word) = size(w.letters)
Base.getindex(w::Word, i::Integer) = w.letters[i]
Base.getindex(w::Word, range::UnitRange) = w.letters[range]
Base.setindex!(w::Word, value, i::Int) = w.letters[i] = value
typeof(::Word{T}) where {T} = T
Base.length(w::Word) = length(w.letters)
getcyclicindex(w::Word, i::Integer) = w.letters[mod1(i, length(w))]
Base.:^(w::AbstractWord, n::Integer) = repeat(w, n)
suffixes(w::AbstractWord) = (w[i:end] for i in firstindex(w):lastindex(w))
function Base.:(==)(letter::T, w::AbstractWord{T}) where {T}
	return length(w) == 1 ? (return w[1] == letter) : false
end
Base.:(==)(w::AbstractWord{T}, letter::T) where {T} = (letter == w)
Base.push!(w::Word{T}, l::T) where {T} = push!(w.letters, l)
Base.first(w::Word) = first(w.letters)
Base.last(w::Word) = last(w.letters)

function isprefix(v::Word, w::Word)
	length(v) ≤ length(w) || return false
	for i ∈ 1:length(v)
		v[i] == w[i] || return false
	end
	return true
end

function issuffix(v::Word, w::Word)
	length(v) <= length(w) || return false
	for i ∈ 1:length(v)
		v[i] == w[i] || return false
	end
	return true
end

"""
    arecyclicallyequal(w::Word{T}, v::Word{T})

	Return whether w is a cyclic permutation of v.
"""
function arecyclicallyequal(w::Word{T}, v::Word{T}) where {T}
	length(w) == length(v) || return false
	for i ∈ 1:length(w)
		outer_continue = false
		for j ∈ 1:length(w)
			if getcyclicindex(w, i - 1 + j) ≠ v[j]
				outer_continue = true
				break
			end
		end
		if outer_continue continue end
		return true
	end
	return false
end

function Base.popfirst!(w::Word)
	return popfirst!(w.letters)
end

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

"""
Replace all occurences of the letter l in w with the word v
"""
function replace_letter!(w::Word{T}, x::T, v::Word{T}) where {T}
	return replace_all!(w, [x], [v])
end

"""
Sequentially replace 
"""
function replace_all!(w::Word{T}, X::Vector{T}, V::Vector{Word{T}}) where {T}
	@assert length(X) == length(V)
	letters = Vector{T}()
	for l ∈ w
		break_outer = false
		for (i, x) ∈ enumerate(X)
			if l == x
				append!(letters, V[i].letters)
				break
			end
		end
	end
	w.letters = letters
	return w
end

function Base.show(io::IO, ::MIME"text/plain", w::AbstractWord)
    if isone(w)
        print(io, 'ε')
    else
        l = length(w)
        for (i, letter) in enumerate(w)
            print(io, letter)
            if i < l
                print(io, '·')
            end
        end
    end
end

function string_repr(w::AbstractWord, A::Alphabet)
    if isone(w)
        return sprint(show, w)  # todo: sprint is probably wrong?
    else
        return join((A[idx] for idx in w), '·')
    end
end

function freeRewriteV1!(w::Word, A::Alphabet)
    i = 1
	@inbounds while i < length(w)
		if hasinverse(A, w[i]) && inv(A, w[i]) == w[i + 1]
			mul!(w, Word(w[begin:(i - 1)]), Word(w[(i + 2):end]))
			i = max(1, i - 1)
		end
		i += 1
	end
	return w
end

function freerewriteBV!(w::Word, A::Alphabet)
	wordlength = length(w)
	wordlength ≥ 2 || return w
	l, r = 1, 2
	mask = trues(wordlength)
	jump = collect(0:wordlength + 1)
	while r ≤ wordlength
		if hasinverse(A, w[l]) && inv(A, w[l]) == w[r]
			mask[l] = false
			mask[r] = false
			jumpnext = jump[l] ≥ 1 ? jump[l] : r += 1
			if l - 1 ≥ 1
				jump[l] = !mask[l - 1] ? jump[l - 1] : l - 1
			end
			jump[r] 	= jump[l]
			jump[r + 1] = jump[l]
			l = jumpnext
		else
			l = r
		end
		r += 1
	end

	letters = w.letters[mask]
	w.letters = letters
	return w
end

cyclically_reduce(w::Word{T}, A::Alphabet{T}) where {T} = cyclically_reduce!(deepcopy(w), A)

function cyclically_reduce!(w::Word{T}, A::Alphabet{T}) where {T}
	freerewriteBV!(w, A)
	i = 1
	while (length(w) - 2 * (i - 1) > 1) && isinverse(A, w[i], w[end - i + 1])
		i += 1
	end
	return Word(w[i:(length(w) - i + 1)])
end

# This is the version from the lecture
function rewrite(
    w::W,
    rewriting,
    vbuffer = one(w),
    wbuffer = one(w),
) where W
	resize!(vbuffer, 0) # empty vbuffer
	
	# copy the content of w to wbuffer, possibly adjusting its size
	resize!(wbuffer, length(w))
	copy!(wbuffer, w)
	
	# do the destructive rewriting from `wbuffer` to `vbuffer`
    v = rewrite!(vbuffer, wbuffer, rewriting)
    return W(v) # return the result of the same type as w
	# and not-aligned with any args passed!
end

function rewrite!(v::Word, w::Word, A::Alphabet)
	for l in w
		if isone(v)
			push!(v, l)
		elseif hasinverse(A, v[end]) && inv(A, v[end]) == l
			resize!(v, length(v) - 1)
		else
			push!(v, l)
		end
	end
    return v
end

macro word_str(string::String)
	length(string) == 0 && return Word{Symbol}()
	letters = [Symbol(s) for s in string]
	return :(Word($letters))
end

struct Words{T}
    A::Alphabet{T}
    wordlength::Int
	iterator

    function Words(A::Alphabet{T}, wordlength::Int) where {T}
		iterator = Iterators.product(ntuple(_ -> A.letters, wordlength)...)
        new{T}(A, wordlength, iterator)
    end
end

function Base.iterate(words::Words)
	iteration = iterate(words.iterator)
	if isnothing(iteration)
		return nothing
	else
		w, initialstate = iteration
		return Word(collect(w)), initialstate
	end
end

function Base.iterate(words::Words, state)
    iteration = iterate(words.iterator, state)
	if isnothing(iteration)
		return nothing
	else
		w, nextstate = iteration
		return Word(collect(w)), nextstate
	end
end

"""
    CyclicWords(A, wordlength)

	Provide an interface to an iterator for some cyclically distinct words over
	given letters in canonical form.

	This is currently not fully implemented. The algorithm misses a lot of words.

	# Examples
	Let A = [a, b] and wordlength = 3. Then the words produced are:
	aaa, aab, abb, bbb.


"""
struct CyclicWords{T}
	letters::Vector{T}
    wordlength::Int
	partition_iterator

	function CyclicWords(A::Alphabet{T}, wordlength::Int) where {T}
		@warn "Doesn't produce all cyclically unique words!"
		partition_iterator = partitions(wordlength, length(A))
		new{T}(A.letters, wordlength, partition_iterator)
	end

	function CyclicWords(A::Vector{T}, wordlength::Int) where {T}
		@warn "Doesn't produce all cyclically unique words!"
		partition_iterator = partitions(wordlength, length(A))
		new{T}(letters, wordlength, partition_iterator)
	end 
end

function word_frompartition(letters::Vector{T},
						    partition::Vector{Int},
							wordlength=sum(partition)::Int) where {T}

	@assert length(letters) == length(partition)
	out = Vector{T}()
	resize!(out, wordlength)
	position = 1
	for (i, l) ∈ enumerate(letters)
		for j ∈ 1:partition[i]
			out[position] = l
			position += 1
		end
	end
	return out
end

function Base.iterate(cyclicwords::CyclicWords)
	partition = iterate(cyclicwords.partition_iterator)
	isnothing(partition) && return nothing
	return iterate(cyclicwords, (partition[1], permutation))
end

function Base.iterate(cyclicwords::CyclicWords, state)
	partition, permutation = state
	if !isnothing(partition)
		if permutation ≤ factorial(big(length(cyclicwords.letters)))
			w = word_frompartition(cyclicwords.letters, partition, cyclicwords.wordlength)
			return Word(w), (partition, permutation + 1)
		else
			nextpartition = iterate(cyclicwords.partition_iterator, partition)
			return iterate(cyclicwords, (nextpartition[1], 1))
		end
	else
		return nothing
	end
end