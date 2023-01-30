import Base.Order: lt, Ordering
include("Alphabet.jl")
include("Word.jl")

abstract type WordOrdering <: Ordering end

struct LenLex{T} <: WordOrdering
    A::Alphabet{T}
    lexicographic_order::AbstractVector{T}
end

function lt(ordering::LenLex, lp::Integer, lq::Integer)
    return ordering.lexicographic_order[lp] < ordering.lexicographic_order[lq]
end

function lt(o::LenLex, p::AbstractWord, q::AbstractWord)
    if length(p) == length(q)
		for (lp, lq) in zip(p, q)
			if lp == lq
				continue
			elseif lt(o, lp, lq)
				return true
			else
				return false
			end
		end
		return false # i.e. p == q
	else
		return length(p) < length(q)
	end
end

