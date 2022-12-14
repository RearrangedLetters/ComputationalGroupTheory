abstract type Group end
abstract type AbstractPermutationGroup{P<:AbstractPermutation} <: Group end

#=
Todo:
    • Introduce dictionary that stores knowledge about the group, such as:
        • isAbelian
        • isSolvable
        • isPolycyclic
        • ...
=#


mutable struct PermutationGroup{P} <: AbstractPermutationGroup{P}
    S::Vector{P}
    order::BigInt
    stabilizerChain::PointStabilizer

    PermutationGroup(S::AbstractVector{P}) where {P<:AbstractPermutation} = new{P}(S)
    PermutationGroup(S::AbstractVector{P}, order::Integer) where {P<:AbstractPermutation} = new{P}(S, order)

    function PermutationGroup(S::AbstractVector{P}, stabilizerChain::PointStabilizer)
        new{P}(S, order(stabilizerChain), schreierSims(S))
    end

    function PermutationGroup(S::AbstractVector{P<:AbstractPermutationGroup}, 
                              order::BigInt, 
                              stabilizerChain::PointStabilizer, 
                              check=true)
        if check
            @assert order == order(stabilizerChain)
            @assert all(S) do s
                _, r = sift(s, stabilizerChain)
                isone(r)
            end
        end
        new{P}(S, order, stabilizerChain)
    end
end

unsafeGenerators(G::PermutationGroup) = G.S
generators(G::PermutationGroup) = copy(G.S)

function stabilizerChain(G::Group)
    if !isdefined(G, :stabilizerChain)
        G.stabilizerChain = schreierSims(G)
        G.order = order(G.stabilizerChain)
    end
    return G.order
end

function order(G::Group)
    if !isdefined(G, :order)
        G.order = order(stabilizerChain(G))
    end
    return G.order
end

function Base.in(G::PermutationGroup, p::AbstractPermutation)
    _, r = sift(getStabilizerChain(G), p)
    return isone(r)
end

Base.one(G::AbstractPermutationGroup) = one(first(getUnsafeGenerators(G)))

Base.eltype(::Type{<:AbstractPermutationGroup{P}}) where P = P

Base.length(G::Group) =
	order(G) > typemax(Int) ? typemax(Int) : order(Int, G)