#=
Todo:
    • Use atomic to make it thread-safe
=#

mutable struct UnionFind
    #=
    Maintains a partition of the set {1:n}
    =#
    S::Set{Any}
    n::Integer
    parents::Vector{UInt}
    ranks::Vector{UInt}
    f

    function UnionFind(n::Integer)
        new(Set(collect(UInt, 1:n)),
            n,
            collect(UInt, 1:n),
            zeros(UInt, n),
            x -> convert(UInt, x))
    end

    function UnionFind(S::Set, f)
        max = maximum([f(s) for s in S])
        new(S, max, collect(UInt, 1:max), zeros(UInt, max), f)
    end
end

Base.length(unionFind::UnionFind) = unionFind.n

function find!(unionFind::UnionFind, i)
    #=
    Returns the representative of i.
    =#
    i = unionFind.f(i)
    if unionFind.parents[i] == i
        return i
    end
    j = find!(unionFind, unionFind.parents[i])
    unionFind.parents[i] = j
    return j
end

function link!(unionFind::UnionFind, i::UInt, j::UInt)
    #=
    Joins two disjoint blocks.
    =#
    @assert find!(unionFind, i) ≠ find!(unionFind, j)
    if unionFind.ranks[i] < unionFind.ranks[j]
        unionFind.parents[i] = j
    else
        unionFind.parents[j] = i
        if unionFind.ranks[i] == unionFind.ranks[j]
            unionFind.ranks[i] += 1
        end
    end
end

function Base.union!(unionFind::UnionFind, i, j)
    #=
    Joins two blocks.
    =#
    rᵢ = find!(unionFind, unionFind.f(i))
    rⱼ = find!(unionFind, unionFind.f(j))
    if rᵢ ≠ rⱼ
        link!(unionFind, rᵢ, rⱼ)
    end
end

function collectBlocks(unionFind::UnionFind)
    blocks = Vector{Vector{eltype(unionFind.S)}}()
    nextBlockIndex::UInt = 1
    blockAssignment = Dict{UInt, Int}()
    for s in unionFind.S
        rᵢ = find!(unionFind, unionFind.f(s))
        if !haskey(blockAssignment, rᵢ)
            blockAssignment[rᵢ] = nextBlockIndex
            nextBlockIndex += 1
        end
        blockIndex = blockAssignment[rᵢ]
        if blockIndex > length(blocks)
            push!(blocks, Vector{eltype(unionFind.S)}())
        end
        push!(blocks[blockIndex], s)
    end
    return blocks
end