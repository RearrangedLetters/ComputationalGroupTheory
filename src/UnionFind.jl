mutable struct UnionFind
    #=
    Maintains a partition of the set {1:n}
    =#
    n::UInt
    parents::Vector{UInt}
    ranks::Vector{UInt}

    function UnionFind(n::Int)
        new(convert(UInt, n), collect(UInt, 1:n), zeros(UInt, n))
    end
end

function find!(unionFind::UnionFind, i::Integer)
    #=
    Returns the representative of i.
    =#
    # i = convert(UInt, i)
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

function Base.union!(unionFind::UnionFind, i::Integer, j::Integer)
    #=
    Joins two blocks.
    =#
    rᵢ = find!(unionFind, convert(UInt, i))
    rⱼ = find!(unionFind, convert(UInt, j))
    if rᵢ ≠ rⱼ
        link!(unionFind, rᵢ, rⱼ)
    end
end

function collectBlocks(unionFind::UnionFind)
    # numberOfBlocks = length(unique(unionFind.parents))
    blocks = Vector{Vector{UInt}}()
    nextBlockIndex::UInt = 1
    blockAssignment = Dict{UInt, Int}()
    for i in 1:(unionFind.n)
        rᵢ = find!(unionFind, i)
        if !haskey(blockAssignment, rᵢ)
            blockAssignment[rᵢ] = nextBlockIndex
            nextBlockIndex += 1
        end
        block = blockAssignment[rᵢ]
        if block > length(blocks)
            push!(blocks, Vector{UInt}())
        end
        push!(blocks[block], i)
    end
    return blocks
end