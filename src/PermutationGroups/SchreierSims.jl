using ComputationalGroupTheory

function Base.:(∈)(g::AbstractPermutation, S::AbstractVector{<:AbstractPermutation})
	return sift(stabilizerChain(S), g) == one(Permutation([1]))
end

function randomlyPickTwo(n::Integer)
    @assert n ≥ 2
    i = rand(1:n)
    j = 0
    while true
        j = rand(1:n)
        j == i || break
    end
    return i, j
end

function pseudorandom(X::AbstractVector{<:GroupElement}, h::GroupElement)
    i, j = randomPick(2, X)
    exponent = random([1, -1])
    side = random([-1, 1])
    if side == -1
        if exponent == -1
            X[i] = inv(X[j]) * X[i]
        else
            X[i] = X[j] * X[i]
        end
        g = X[i] * h
    else
        if exponent == -1
            X[i] = X[i] * inv(X[j])
        else
            X[i] = X[i] * X[j]
        end
        g = h * X[i]
    end
    return X ∪ [g], g
end

@inline function pad(S::AbstractVector{<:GroupElement}, length::Integer)
    Δ = length - length(S)
    return Δ > 0 ? append!(S, fill(Permutation([1]), Δ)) : S
end

function pseudorandomList(S::AbstractVector{<:GroupElement}, n::Integer=50)
    X = pad(S, 11)
    h = Permutation([1])
    for _ in 1:n
        X, h = pseudorandom(X, h)
    end
    return X
end

mutable struct PointStabilizer{P<:AbstractPermutation}
    S::AbstractVector{P}
    T::AbstractTransversal
    stabilizer::PointStabilizer{P}

    PointStabilizer{P}() where P = new{P}(Vector{P}())
    #= function PointStabilizer{P}(transversal::AbstractTransversal) where P
        new{P}(Vector{P}(), transversal)
    end =#

    function PointStabilizer{P}(stabilizer::PointStabilizer) where P
        new{P}(stabilizer.S, stabilizer.T, stabilizer.stabilizer)
    end
end

function copy(stabilizer::PointStabilizer{P}) where P
    return PointStabilizer{P}(stabilizer)
end

#= mutable struct PointStabilizer{P<:AbstractPermutation, R<:AbstractOrbit}
    S::AbstractVector{P}
    T::R
    stabilizer::PointStabilizer{P, R}

    function PointStabilizer{P::AbstractPermutation, R::AbstractOrbit}()
        PointStabilizer{P, R}() where P = new{P, R}(Vector{P}())
    end
end =#

generators(pointStabilizer::PointStabilizer) = pointStabilizer.S
point(pointStabilizer::PointStabilizer) = pointStabilizer.T.x
transversal(pointStabilizer::PointStabilizer) = pointStabilizer.T
stabilizer(pointStabilizer::PointStabilizer) = pointStabilizer.stabilizer
Base.isempty(pointStabilizer::PointStabilizer) = isempty(generators(pointStabilizer))

function schreierSims(S::AbstractVector{<:AbstractPermutation})
    @assert !isempty(S)
    pointStabilizer = PointStabilizer{eltype(S)}()
    for s in S
        push!(pointStabilizer, s)
    end
    return pointStabilizer
end

function Base.push!(pointStabilizer::PointStabilizer, g::AbstractPermutation)
    g = sift(pointStabilizer, g)
    if isone(g)
        return pointStabilizer
    end
    if isempty(pointStabilizer)
        extendChain!(pointStabilizer, g)
    else
        extendGenerators!(pointStabilizer, g)
    end

    return pointStabilizer
end

function sift(pointStabilizer::PointStabilizer, g::AbstractPermutation, L = AbstractPermutation[])
    # returns 1 iff g is in the point stabilizer
    if isempty(pointStabilizer) || isone(g)
        return g, L
    else
        x = point(pointStabilizer)
        δ = x^g
        T = transversal(pointStabilizer)
        if δ ∈ T
            g = g * inv(T[δ])  # point in the stabilizer of x
            @assert x^g == x
            return sift(stabilizer(pointStabilizer), g, push!(T[δ]))
        else
            return g, L
        end
    end
end

@inline firstMoved(g::AbstractPermutation) = findfirst(x -> x^g ≠ x, 1:(degree(g) + 1))

function order(g::AbstractPermutation)
    n = 1
    e = one(g)
    h = deepcopy(g)
    while h^n ≠ e
        h *= g
        n += 1
    end
    return n
end

function extendChain!(pointStabilizer::PointStabilizer{P}, g::AbstractPermutation) where P
    @assert !isone(g)
    push!(pointStabilizer.S, g)
    pointStabilizer.T = FactoredTransversal(generators(pointStabilizer), firstMoved(g))  # todo: specify type?
    pointStabilizer.stabilizer = PointStabilizer{P}()

    k = length(transversal(pointStabilizer))
    if k < order(g)
        extendChain!(stabilizer(pointStabilizer), g^k)
    end

    return pointStabilizer
end

function extendGenerators!(pointStabilizer::PointStabilizer, g::AbstractPermutation)
    @assert !isone(g)
    # simple version
    push!(pointStabilizer.S, g)
    pointStabilizer.T = FactoredTransversal(generators(pointStabilizer), point(pointStabilizer))
    T = transversal(pointStabilizer)

    for s in generators(pointStabilizer)
        for δ in transversal(pointStabilizer)  # iteration over points in the orbit
            r = T[δ]
            schreierGenerator = r * s * inv(T[δ^s])
            if !isone(schreierGenerator)
                push!(stabilizer(pointStabilizer), schreierGenerator)
            end
        end
    end
    return pointStabilizer
end

function Base.iterate(pointStabilizer::PointStabilizer,
                      iterator::PointStabilizer=pointStabilizer)

    if isdefined(iterator, :stabilizer)
        return iterator, iterator.stabilizer
    else
        return nothing
    end
end

function order(S::AbstractVector{<:GroupElement})
    𝒞 = schreierSims(S)
    order = BigInt(1)
    for C in 𝒞
        order *= length(transversal(C))
    end
    return order
end

function Base.show(io::IO, pointStabilizer::PointStabilizer)
    if isempty(pointStabilizer)
        print(io, "⟨⟩")
    else
        for C in pointStabilizer
            if isempty(C)
                print(io, "⟨⟩")
            else
                print(io, "G_i = ⟨$(generators(C))⟩, ")
                print(io, "x = $(point(C)), ")
                print(io, "T = $(transversal(C))")
            end
        end
    end
end

function Base.length(pointStabilizer::PointStabilizer)
    count = 0
    for _ in pointStabilizer
        count += 1
    end
    return count
end

function isAbelian(S::AbstractVector{Permutation})
    for s in S
        for t in S
            if s^t ≠ t^s
                return false
            end
        end
    end
    return true
end

#=
struct StabilizerChain
    S::AbstractVector{AbstractVector{<:GroupElement}}
    β::AbstractVector{Integer}
    T::Transversal
end

function push!(stabilizerChain::StabilizerChain, g::AbstractPermutation, d::Integer)
    for i in 1:(d - 1)
        @assert stabilizerChain.β[i]^g == stabilizerChain.β[i]
    end

    if d > length(stabilizerChain)
        β = firstMoved(g)
        S = [g]
        T = Transversal(β, S)
        extendBase(stabilizerChain, β)
        extendTransversal(stabilizerChain, T)
        k = length(T)
        if k < order(g)
            push!(stabilizerChain, g^k, d + 1)  # todo: g^k is not "correct"
        end
    else
        push!(stabilizerChain.S[d], g)
        for s in schreier(stabilizerChain.β[d], stabilizerChain.S[d])  # todo: correct call of schreier?
            _, r = sift(stabilizerChain, s)
            if r ≠ one(first(S))
                push!(stabilizerChain, s, d + 1)
            end
        end
    end
end


=#