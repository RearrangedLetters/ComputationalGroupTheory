abstract type Group end
abstract type AbstractPermutationGroup{P<:AbstractPermutation} <: Group end

mutable struct PermutationGroup{P} <: AbstractPermutationGroup{P}
    S::Vector{P}
    order::BigInt
    stabilizerChain::PointStabilizer
    knowledge::Dict{String, Union{Bool, Nothing}}(
        "order"         => false,
        "abelian"       => nothing,
        "solveable"     => nothing,
        "polycyclic"    => nothing,)

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
        G.stabilizerChain = if knowsOrder(G) schreierSims(G, order(G)) else schreierSims(G) end
    end
    return G.stabilizerChain
end

knowsOrder(G::Group) = !(G.knowledge["order"] === nothing)

function order(G::Group)
    if !knowsOrder(G)
        G.order = order(stabilizerChain(G))
        G.knowledge["order"] = true
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
	order(G) > typemax(Int) ? typemax(Int) : convert(Int, order(G))

function isAbelian(G::AbstractPermutationGroup)
    if G.knowledge["abelian"] === nothing
        for s in G.S
            for t in G.S
                if s^t ≠ t^s
                    G.knowledge["abelian"] = false
                end
            end
        end
        G.knowledge["abelian"] = true
    end
    return G.knowledge["abelian"]
end