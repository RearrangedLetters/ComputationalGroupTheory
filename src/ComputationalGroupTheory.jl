module ComputationalGroupTheory

export Permutation, degree, orbit, transversal, transversalFactored, transversalSchreier, representative,
       @perm_str, schreierSims, order, PointStabilizer, Transversal, FactoredTransversal,
       TransversalTree, enumerateGroup, UnionFind, find, union, collectBlocks, nullspace,
       echelonize

include("AbstractPermutation.jl")
include("CyclePermutation.jl")
include("Permutation.jl")
include("OrbitTransversalSchreier.jl")
include("SchreierSims.jl")
include("Transversal.jl")
# include("Backtrack.jl")
include("UnionFind.jl")
include("ExactMatrices.jl")

end