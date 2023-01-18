module ComputationalGroupTheory

using Random

export Permutation, degree, orbit, transversal, transversalFactored, transversalSchreier, representative,
       @perm_str, schreierSims, order, PointStabilizer, Transversal, FactoredTransversal,
       TransversalTree, backtrackRecursive, backtrack, isAbelian

       
include("PermutationGroups/AbstractPermutation.jl")
include("PermutationGroups/CyclePermutation.jl")
include("PermutationGroups/Permutation.jl")
include("PermutationGroups/OrbitTransversalSchreier.jl")
include("PermutationGroups/SchreierSims.jl")
include("PermutationGroups/Transversal.jl")
include("PermutationGroups/Backtrack.jl")

include("Algorithms/UnionFind.jl")
export UnionFind, find, union, collectBlocks

include("Cohomology/ExactMatrices.jl")
export nullspace, rowspace, echelonize, inv

# include("Group.jl")

include("FreeGroups/Alphabet.jl")
export Alphabet, setinverse!, hasinverse

include("Experimenting/Residue.jl")
export Residue, zero

include("Cohomology/Cohomology.jl")
export Z1, Z1!, B1, B1!

include("FreeGroups/Word.jl")
export AbstractWord, Word, one, freeRewriteBV!, rewrite, @Σ_str, @w_str

end