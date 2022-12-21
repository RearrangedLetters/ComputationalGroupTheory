import Primes

struct Residue
    residue::Int
    modulus::Int

    function Residue(residue, modulus)
        @assert modulus > 1 "Modulus must be greater than 1, got $modulus"
        return new(mod(residue, modulus), modulus)
    end
end

Base.show(io::IO, x::Residue) = print(io, x.residue, " mod ", x.modulus)

function Base.:+(x::Residue, y::Residue)
    @assert x.modulus == y.modulus "Addition of residues with different modulos is not well-
                                    defined!"
    return Residue(x.residue + y.residue, x.modulus)
end

function Base.:-(x::Residue)
    return Residue(-x.residue, x.modulus)
end

function Base.:-(x::Residue, y::Residue)
    @assert x.modulus == y.modulus "Addition of residues with different modulos is not well-
                                    defined!"
    return x + (-y)
end

function Base.zero(x::Residue)
    return Residue(0, x.modulus)
end

function Base.:(==)(x::Residue, y::Residue)
    return x.residue == y.residue && x.modulus == y.modulus
end

function Base.:+(x::Residue, y::Int)
    return Residue(x.residue + y, x.modulus)
end

function Base.:+(x::Int, y::Residue)
    return Residue(x + y.residue, y.modulus)
end

function Base.:+(x::BigInt, y::Residue)
    return Residue(x + y.residue, y.modulus)
end

function Base.:*(x::Residue, y::Residue)
    @assert x.modulus == y.modulus "Multiplication of residues with different modulos is not
                                    well-defined!"
    return Residue(x.residue * y.residue, x.modulus)
end

function Base.inv(x::Residue)
    return Residue(x.residue^(Primes.totient(x.modulus) - 1), x.modulus)
end