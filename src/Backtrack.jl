include("SchreierSims.jl")

"""
    Vertex{V}
Abstract type that defines a single function to get the content of type V of a vertex:
    • get(v::Vertex)
"""
abstract type AbstractVertex{V} end

"""
    AbstractTree{L}
Abstract type with interface to work with trees. Vertices are of type V.
The only method of the interface is
    • getRoot()
"""
abstract type AbstractTree{V} end

"""
    AbstractBacktrackTree
Defines an abstract type to work with implicit trees. I.e. the tree isn't held in memory,
but instead only is locally procuded by the methods of the interface. Besides getRoot() from the
super-type, there are two more methods to implement:
    • getChild(v::Vertex)
    • nextSibling(v::Vertex)
    • isLeaf(v::Vertex)
"""
abstract type AbstractBacktrackTree{V} <: AbstractTree{V} end

struct TransversalTree{V} <: AbstractBacktrackTree{V}
    𝒞::PointStabilizer{}

    function TransversalTree(S::AbstractVector{<:AbstractPermutation}) where V
        𝒞 = schreierSims(S)
        new{eltype(first(𝒞.T))}(𝒞)
    end
end

function transversal(transversalTree::TransversalTree, depth::Int)
    @assert depth ≤ length(transversalTree.𝒞)
    i::Int = 1
    pointStabilizer::PointStabilizer = transversalTree.𝒞
    while i < depth
        pointStabilizer = pointStabilizer.stabilizer
        i += 1
    end
    return pointStabilizer.T
end

function backtrackRecursive(S::AbstractVector{<:AbstractPermutation}, oracle=(x -> true))
    return backtrackHelper!(Vector{eltype(S)}(), 
                                 oracle,
                                 TransversalTree(S),
                                 one(first(S)),
                                 1)
end

function backtrackHelper!(L::AbstractVector{<:AbstractPermutation},
                          oracle,
                          tree::TransversalTree,
                          g::AbstractPermutation,
                          depth::Int)
    if oracle(g)
        T = transversal(tree, depth)
        for δ in T
            if length(tree.𝒞) == depth
                push!(L, g * T[δ])
            else
                backtrackHelper!(L, oracle, tree, g * T[δ], depth + 1)
            end
        end
    end
    return L
end

function getStabilizer(stabilizer::PointStabilizer, depth::Int)
    i::Int = 1
    pointStabilizer = copy(stabilizer)
    while i < depth
        pointStabilizer = pointStabilizer.stabilizer
        i += 1
    end
    return pointStabilizer
end

function getTransversal(stabilizer::PointStabilizer, depth::Int)
    return getStabilizer(stabilizer, depth).T
end

function depth(stabilizer::PointStabilizer)
    d = 0
    while isdefined(stabilizer, :stabilizer)
        stabilizer = stabilizer.stabilizer
        d +=1
    end
    return d
end

function last(stabilizer::PointStabilizer)
    stabilizer = copy(stabilizer)
    while isdefined(stabilizer, :stabilizer)
        stabilizer = stabilizer.stabilizer
    end
    return stabilizer
end

function backtrack(S::AbstractVector{<:AbstractPermutation})
    stabilizerChain = schreierSims(S)
    s = Vector{eltype(stabilizerChain.T.Ωᴳ)}()
    indices = Int[1]
    d = depth(stabilizerChain)
    g = [one(first(S))]
    while true
        stabilizer = getStabilizer(stabilizerChain, length(indices))
        T = getTransversal(stabilizerChain, length(indices))
        if Base.last(indices) <= length(T)
            δ = T.Ωᴳ[Base.last(indices)]
            push!(s, δ^Base.last(g))
            if length(indices) < d
                push!(g, T[δ])
                push!(indices, 1)
            else
                for δ ∈ T
                    println(T[δ] * foldl(*, g))
                end
                pop!(s), pop!(s), pop!(indices), pop!(g)
                indices[length(indices)] += 1
            end
        else
            break
        end
    end
end

function Base.iterate(S::AbstractVector{<:AbstractPermutation}, 
                      i=([], [], [1], [one(first(S))]))
    stabilizerChain = if isempty(i[1]) schreierSims(S) else i[1] end
    s = i[2]
    indices = i[3]
    d = depth(stabilizerChain)
    g = i[4]
    while true
        stabilizer = getStabilizer(stabilizerChain, length(indices))
        T = getTransversal(stabilizerChain, length(indices))
        if Base.last(indices) <= length(T)
            δ = T.Ωᴳ[Base.last(indices)]
            push!(s, δ^Base.last(g))
            if length(indices) < d
                push!(g, T[δ])
                push!(indices, 1)
            else
                for δ ∈ T
                    indices[length(indices) - 1] += 1
                    return (T[δ] * foldl(*, g), 
                           (stabilizerChain, 
                            s[:length(s) - 2],
                            indices[:length(indices) - 1], g[:length(g) - 1]))
                end
            end
        else
            return nothing
        end
    end
end