import Base.Order: lt, Ordering

abstract type WordOrdering <: Ordering end

struct LenLex{T} <: WordOrdering
    alphabet::Alphabet{T}
    lexicographic_ordering::AbstractVector{T}
end

function lt(ordering::LenLex, lp::Integer, lq::Integer)
    return ordering.lexicographic_ordering[lp] < ordering.lexicographic_ordering[lq]
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

