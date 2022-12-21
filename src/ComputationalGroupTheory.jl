module ComputationalGroupTheory

using Random

export Permutation, degree, orbit, transversal, transversalFactored, transversalSchreier, representative,
       @perm_str, schreierSims, order, PointStabilizer, Transversal, FactoredTransversal,
       TransversalTree, backtrackRecursive, backtrack, UnionFind, find, union, collectBlocks, nullspace,
       echelonize, isAbelian

export Alphabet, setinverse!, hasinverse

include("AbstractPermutation.jl")
include("CyclePermutation.jl")
include("Permutation.jl")
include("OrbitTransversalSchreier.jl")
include("SchreierSims.jl")
include("Transversal.jl")
include("Backtrack.jl")
include("UnionFind.jl")
include("ExactMatrices.jl")
# include("Group.jl")
include("Alphabets.jl")

end