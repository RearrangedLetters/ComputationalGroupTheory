include("../PermutationGroups/Group.jl")

mutable struct MatrixGroup{T} <: Group
    S::Vector{Matrix{T}}

    function MatrixGroup(S::Vector{Matrix{T}}) where {T}
        new{T}(S)
    end
end

Base.one(G::MatrixGroup{T}) where {T} = one(Matrix{T})