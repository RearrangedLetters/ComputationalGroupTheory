abstract type GroupElement end

"""
AbstractPermutation
Abstract type representing permutations of set `1:n`.

Subtypes `Perm <: AbstractPermutation` must implement the following functions:
* `(σ::Perm)(i::Integer)` - the image of `i` under `σ`,
* `degree(σ::Perm)` - the minimal `n` such that `σ(k) == k` for all `k > n`,
* `Perm(images::AbstractVector{<:Integer}[, check::Bool=true])` - construct a
`Perm` from a vector of images. Optionally the second argument `check` may be
set to `false` when the caller knows that `images` constitute an honest
permutation.
"""
abstract type AbstractPermutation <: GroupElement end

function degree end

Base.one(σ::P) where P<:AbstractPermutation = P(Int[], false)
Base.isone(σ::AbstractPermutation) = degree(σ) == 1

function Base.inv(σ::P) where P<:AbstractPermutation
    images = similar(1:degree(σ))
    for i in 1:degree(σ)
        images[i^σ] = i
    end
    return P(images, false)
end

function Base.:(*)(σ::P, τ::AbstractPermutation) where P<:AbstractPermutation
    aDegree = max(degree(σ), degree(τ))
    images = similar(1:aDegree)
    for i in 1:aDegree
        images[i] = (i^σ)^τ
    end
    return P(images, false)
end

function Base.:(==)(σ::AbstractPermutation, τ::AbstractPermutation)
    degree(σ) ≠ degree(τ) && return false
    for i in 1:degree(σ)
        if i^σ != i^τ
            return false
        end
    end
    return true
end

function Base.hash(σ::AbstractPermutation, h::UInt)
    h = hash(AbstractPermutation, h)
    for i in 1:degree(σ)
        h = hash(i^σ, h)
    end
    return h
end
