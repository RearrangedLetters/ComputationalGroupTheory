module ComputationalGroupTheory

using Random

export Permutation, degree, orbit, transversal, transversalFactored, transversalSchreier, representative,
       @perm_str, schreierSims, order, PointStabilizer, Transversal, FactoredTransversal,
       TransversalTree, backtrackRecursive, backtrack, UnionFind, find, union, collectBlocks, nullspace,
       echelonize, isAbelian

       
include("PermutationGroups/AbstractPermutation.jl")
include("PermutationGroups/CyclePermutation.jl")
include("PermutationGroups/Permutation.jl")
include("PermutationGroups/OrbitTransversalSchreier.jl")
include("PermutationGroups/SchreierSims.jl")
include("PermutationGroups/Transversal.jl")
include("PermutationGroups/Backtrack.jl")
include("Algorithms/UnionFind.jl")
include("Cohomology/ExactLinearAlgebra.jl")
# include("Group.jl")

include("FreeGroups/Alphabets.jl")
export Alphabet, setinverse!, hasinverse

end