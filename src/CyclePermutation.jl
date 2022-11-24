include("AbstractPermutation.jl")

struct CyclePermutation <: AbstractPermutation
    cycles::Vector{Vector{Int}}

    function CyclePermutation(images::Vector{<:Integer}, check=true)
        new(cycle_decomposition(Permutation(images, check)))
    end

    function degree(σ::CyclePermutation)
		degree = 1
		for cycle in σ.cycles
			degree = max(degree, max(cycle))
		end
		return degree
	end

    function apply(cycle::Vector{Int}, n::Int)
		for i in 1:length(cycle)
			if cycle[i] == n
				return cycle[i % length(cycle) + 1]
			end
		end
		return n
	end

    function Base.:^(i::Integer, σ::CyclePermutation)
        n = i
		for j in degree(σ):-1:1
			n = apply(σ.cycles[j], n)
		end
		return n
	end
end