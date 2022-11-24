function Base.:(∈)(g::AbstractPermutation, S::AbstractVector{<:AbstractPermutation})
	return sift(stabilizerChain(S), g) == one(Permutation([1]))
end

function normalClosure(S::AbstractVector{<:GroupElement}, U::AbstractVector{<:GroupElement})
    N = copy(U)
    for n in N
        for s in S
            γ = inv(s) * n * s
            if γ ∉ N  # todo: this works if Julia automatically uses ∈ to define ∉
                push!(N, γ)
            end
        end
    end
    return N
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
    x::Integer
    T#::Transversal
    stabilizer::PointStabilizer{P}

    PointStabilizer{P}() where P = new{P}(Vector{P}())  # incomplete initialization

    generators(pointStabilizer::PointStabilizer) = pointStabilizer.S
    # point(pointStabilizer::PointStabilizer) = pointStabilizer.x
    point(pointStabilizer::PointStabilizer) = first(transversal(pointStabilizer))
    transversal(pointStabilizer::PointStabilizer) = pointStabilizer.transversal
    stabilizer(pointStabilizer::PointStabilizer) = pointStabilizer.stabilizer

    Base.isempty(pointStabilizer::PointStabilizer) = isempty(generators(pointStabilizer))  # or: point(pointStabilizer) == 0
end

struct Transversal
    x::Integer
    T::AbstractVector{Integer}

    function Transversal(x::Integer, S::AbstractVector{<:GroupElement})
        _, aTransversal = transversal(S, [x])
        new(x, aTransversal)
    end

    point(transversal::Transversal) = transversal.x
    Base.getindex(transversal::Transversal, i) = transversal.T[i]
    Base.length(transversal::Transversal) = length(transversal.T)
end

function stabilizerChain(S::AbstractVector{<:AbstractPermutation})
    𝒞 = PointStabilizer{eltype(S)}()
    for s ∈ S
        _, r = sift(𝒞, s)
        if r ≠ one(first(S))
            push!(𝒞, r)
        end
    end
    return 𝒞
end

function schreierSims(S::AbstractVector{<:AbstractPermutation})
    @assert !isempty(S)
    pointStabilizer = PointStabilizer{eltype(S)}()
    for s in S
        push!(pointStabilizer, s)
    end
    return pointStabilizer
end

function push!(pointStabilizer::PointStabilizer, g::AbstractPermutation)
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

function sift(pointStabilizer::PointStabilizer, g::AbstractPermutation)  # todo: should return more
    # returns 1 iff g is in the point stabilizer
    if isempty(pointStabilizer) || isone(g)
        return g
    else
        x = point(pointStabilizer)
        δ = x^g
        T = transversal(pointStabilizer)
        if δ in T
            r = T[δ]
            g = g * inv(r)  # point in the stabilizer of x
            @assert x^g == x
            return sift(stabilizer(pointStabilizer), g)
        else
            return g
        end
    end
end

@inline firstMoved(g::AbstractPermutation) = findfirst(x -> x^g ≠ x, 1:(degree(g) + 1))

function extendChain!(pointStabilizer::PointStabilizer{P}, g::AbstractPermutation) where P
    @assert !isone(g)
    push!(pointStabilizer.S, g)
    pointStabilizer.T = Transversal(firstMoved(g), generators(pointStabilizer))
    pointStabilizer.stabilizer = PointStabilizer{P}()

    k = length(pointStabilizer.T)
    if k < order(g)
        extendChain!(stabilizer(pointStabilizer), g^k)
    end

    return pointStabilizer
end

function extendGenerators(pointStabilizer::PointStabilizer, g::AbstractPermutation)
    @assert !isone(g)
    # simple version
    push!(pointStabilizer.S, g)
    pointStabilizer.T = Transversal(point(pointStabilizer), generators(pointStabilizer))
    T = transversal(pointStabilizer)
    for s in generators(pointStabilizer)
        for δ in transversal(pointStabilizer)  # iteration over points in the orbit
            r = T[δ]
            schreier_generator = r * s * inv(T[δ^s])
            if !isone(schreier_generator)
                push!(stabilizer(pointStabilizer), schreier_generator)
            end
        end
    end
    return pointStabilizer
end