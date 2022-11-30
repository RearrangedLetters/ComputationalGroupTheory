"""
    AbstractOrbit{X}
Abstract type representing abstract orbits of elements of type `X`.
"""
abstract type AbstractOrbit{X} end

"""
    AbstractTransversal{X, Y} <: AbstractOrbit{X}
Abstract type representing the bijection of orbit oand orbit representatives.

`X` is the type of elements in the orbit, while `Y` is the type of the
representatives. When `tr` is a transversal of `x` and `g` is a `GroupElement`
then `tr[x^g]` returns the representative of the `g`-coset of the stabilizer of `x`.

## Methods to implement:
 * Constructors:
  - `Transversal(x, g::GroupElement[, action=^])` a specific constructor for a
    cyclic group
  - `Transversal(x, S::AbstractVector{<:GroupElement}[, action=^])` the default
    constructor
 * `Base.getindex(tr::T, n::Integer)` - return the coset representative of
   corresponding to `n`, i.e. a group element `g` such that `first(tr)^g == n`.
   If no such element exists a `NotInOrbit` exception will be thrown.
 * Iteration protocol, iterating over points in the orbit.
"""
abstract type AbstractTransversal{X<:GroupElement, Y} <: AbstractOrbit{X} end

Base.eltype(::Type{<:AbstractTransversal{X}}) where X = X

struct NotInOrbit <: Exception
    x
    first
end

function Base.showerror(io::IO, error::NotInOrbit)
    print(io, error.x, " was not found in the orbit of ", error.first)
end

struct Transversal{X, Y} <: AbstractTransversal{X, Y}
    x::Y
    Ωᴳ::AbstractVector{Y}
    T::AbstractDict{Y, X}

    function Transversal(g::GroupElement, x, action=^)
        Ωᴳ, T = transversal([g], [x], action)
        new{GroupElement, typeof(x)}(x, Ωᴳ, T)
    end

    function Transversal(S::AbstractVector{X}, x, action=^) where X<:GroupElement
        Ωᴳ, T = transversal(S, x, action)
        new{X, typeof(x)}(x, Ωᴳ, T)
    end
end

function Base.getindex(transversal::Transversal, n::Integer)
    if haskey(transversal.T, n)
        return transversal.T[n]
    else
        throw(NotInOrbit(n, transversal.Ωᴳ))
    end
end

Base.length(transversal::Transversal) = length(transversal.Ωᴳ)

function Base.iterate(transversal::Transversal, i=1)
    return i > length(transversal) ? nothing : (transversal.Ωᴳ[i], i + 1)
end

function Base.show(io::IO, transversal::Transversal)
    println(io, transversal.T)
end

struct FactoredTransversal{X, Y} <: AbstractTransversal{X, Y}
    x::Y
    Ωᴳ::AbstractVector{Y}
    T::AbstractDict{Y, <:AbstractVector{X}}

    function FactoredTransversal(g::GroupElement, x, action=^)
        Ωᴳ, T = transversalFactored([g], [x], action)
        new{GroupElement, typeof(x)}(x, Ωᴳ, T)
    end

    function FactoredTransversal(S::AbstractVector{X}, x, action=^) where X<:GroupElement
        Ωᴳ, T = transversalFactored(S, x, action)
        new{X, typeof(x)}(x, Ωᴳ, T)
    end
end

function Base.getindex(transversal::FactoredTransversal, n::Integer)
    if haskey(transversal.T, n)
        g = foldl(*, transversal.T[n])
        transversal.T[n] = [g]
        return g
    else
        throw(NotInOrbit(n, transversal.Ωᴳ))
    end
end

Base.length(transversal::FactoredTransversal) = length(transversal.Ωᴳ)

function Base.iterate(transversal::FactoredTransversal, i=1)
    return i > length(transversal) ? nothing : (transversal.Ωᴳ[i], i + 1)
end