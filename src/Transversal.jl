"""
    AbstractOrbit{S}
Abstract type representing abstract orbits of elements of type `S`.
"""
abstract type AbstractOrbit{S} end

"""
    AbstractTransversal{S, T} <: AbstractOrbit{S}
Abstract type representing the bijection of orbit oand orbit representatives.

`S` is the type of elements in the orbit, while `T` is the type of the
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
abstract type AbstractTransversal{S, T<:GroupElement} <: AbstractOrbit{S} end

Base.eltype(::Type{<:AbstractTransversal{S}}) where S = S

struct NotInOrbit <: Exception
    x
    first
end

function Base.showerror(io::IO, e::NotInOrbit)
    print(io, e.x, " was not found in the orbit of ", e.first)
end

struct Transversal{S, T} <: AbstractTransversal{S, T}
    Ωᴳ::AbstractVector{S}
    T::AbstractDict{S, T}

    function Transversal(g::GroupElement, x, action=^)
        Ωᴳ, T = transversal([g], [x], action)
        typeofS = typeof(x)
        typeofT = typeof(last(first(T)))
        new{typeofS, typeofT}(Ωᴳ, T)
    end

    function Transversal(S::AbstractVector{<:GroupElement}, x, action=^)
        Ωᴳ, T = transversal(S, x, action)
        typeofS = typeof(x)
        typeofT = typeof(last(first(T)))
        new{typeofS, typeofT}(Ωᴳ, T)
    end
end

function Base.getindex(transversal::Transversal, n::Integer)
    if haskey(transversal.T, n)
        return transversal.T[n]
    else
        throw(NotInOrbit(n, n))  # todo: implement this according to the interface
    end
end

Base.length(transversal::Transversal) = length(transversal.Ωᴳ)

function Base.iterate(transversal::Transversal, i=1)
    return i ≥ length(transversal) ? nothing : (transversal[i], i + 1)
end