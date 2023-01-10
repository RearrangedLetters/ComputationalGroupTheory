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

include("FreeGroups/Alphabets.jl")
export Alphabet, setinverse!, hasinverse

include("Experimenting/Residue.jl")
export Residue, zero

include("Cohomology/Cohomology.jl")
export Z1, Z1!, B1, B1!

include("Algorithms/MachineLearning/CNN.jl")
export DenseLayer, Tanh, mean_square_loss, mean_square_loss_prime, forward, backward

end