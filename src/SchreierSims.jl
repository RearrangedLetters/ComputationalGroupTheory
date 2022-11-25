include("Transversal.jl")

function Base.:(∈)(g::AbstractPermutation, S::AbstractVector{<:AbstractPermutation})
	return sift(stabilizerChain(S), g) == one(Permutation([1]))
end

function normalClosure(S::AbstractVector{<:GroupElement}, U::AbstractVector{<:GroupElement})
    N = copy(U)
    for n in N
        for s in S
            γ = inv(s) * n * s
            if γ ∉ N
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
    T::Transversal
    stabilizer::PointStabilizer{P}

    PointStabilizer{P}() where P = new{P}(Vector{P}())  # incomplete initialization
end

generators(pointStabilizer::PointStabilizer) = pointStabilizer.S
point(pointStabilizer::PointStabilizer) = first(transversal(pointStabilizer))
transversal(pointStabilizer::PointStabilizer) = pointStabilizer.T
stabilizer(pointStabilizer::PointStabilizer) = pointStabilizer.stabilizer
Base.isempty(pointStabilizer::PointStabilizer) = isempty(generators(pointStabilizer))

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

function sift(pointStabilizer::PointStabilizer, g::AbstractPermutation)  # todo: should return more according to lecture (but not according to notebook)
    # returns 1 iff g is in the point stabilizer
    if isempty(pointStabilizer) || isone(g)
        return g
    else
        x = point(pointStabilizer)
        δ = x^g
        T = transversal(pointStabilizer)
        if δ in T
            g = g * inv(T[δ])  # point in the stabilizer of x
            @assert x^g == x
            return sift(stabilizer(pointStabilizer), x)
        else
            return g
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

Base.:(^)

function extendChain!(pointStabilizer::PointStabilizer{P}, g::AbstractPermutation) where P
    @assert !isone(g)
    Base.push!(pointStabilizer.S, g)
    pointStabilizer.T = Transversal(generators(pointStabilizer), firstMoved(g))
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
    pointStabilizer.T = Transversal(generators(pointStabilizer), point(pointStabilizer))
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