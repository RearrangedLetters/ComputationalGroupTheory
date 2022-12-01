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
    # T::AbstractTransversal
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

function enumerateGroup(S::AbstractVector{<:AbstractPermutation})
    return enumerateGroupHelper!(Vector{eltype(S)}(), TransversalTree(S), one(first(S)), 1)
end

function enumerateGroupHelper!(L::AbstractVector{<:AbstractPermutation}, tree::TransversalTree,
                          g::AbstractPermutation, depth::Int)
    T = transversal(tree, depth)
    for δ in T
        if length(tree.𝒞) == depth
            push!(L, g * T[δ])
        else
            enumerateGroupHelper!(L, tree, g * T[δ], depth + 1)
        end
    end
    return L
end