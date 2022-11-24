module ComputationalGroupTheory

export Permutation, degree, orbit
include("AbstractPermutation.jl")
include("CyclePermutation.jl")
include("Permutation.jl")
include("OrbitTransversalSchreier.jl")

end