abstract type GroupElement end
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

function Base.show(io::IO, σ::AbstractPermutation)
    if isone(σ)
        print(io, "()")
    else
        for cycle in cycle_decomposition(σ)
            if length(cycle) == 1
                continue
            else
                print(io, "(")
                join(io, cycle, ",")
                print(io, ")")
            end
        end
    end
end

function cycle_decomposition(σ::AbstractPermutation)
    visited = falses(degree(σ))
    cycles = Vector{Vector{Int}}()
    # each cycle will be a Vector{Int} and we have a whole bunch of them
    for i in 1:degree(σ)
        if visited[i]
            # if we have already seen this point there is no point in computing
            # the same orbit twice
            continue # i.e. skip the rest of the body and continue with the next i
        end
        Δ = orbit_plain(i, σ, ^)
        visited[Δ] .= true # modify the `visited` along the whole orbit
        push!(cycles, Δ) # add obtained orbit to cycles
    end
    return cycles
end

struct Permutation <: AbstractPermutation
    images::Vector{UInt}

    function Permutation(v::Vector{<:Integer}, check=true)
        if check
            @assert sort(v) == 1:length(v) "Image vector doesn't define a permutation."
        end
        return new(v)
    end

    function Base.:^(i::Integer, σ::Permutation)
        return σ.images[i]
    end

    function degree(σ::Permutation)
        return something(findlast(i -> σ.images[i] != i, 1:length(σ.images)), 1)
    end
end

struct CyclePermutation <: AbstractPermutation
    cycles::Vector{Vector{Int}}

    function CyclePermutation(images::Vector{<:Integer}, check=true)
        new(cycle_decomposition(Permutation(images, check)))
    end

    function degree(σ::CyclePermutation)
		degree = 1
		for cycle in σ.cycles
			degree = max(degree, max(cycle))
		end
		return degree
	end

    function apply(cycle::Vector{Int}, n::Int)
		for i in 1:length(cycle)
			if cycle[i] == n
				return cycle[i % length(cycle) + 1]
			end
		end
		return n
	end

    function Base.:^(i::Integer, σ::Permutation)
		for i in length(σ):-1:1
			n = apply(σ.cycles[i], n)
		end
		return n
	end
end

function orbit(s::GroupElement, ω, action=^)
    #=
    In: • G = ⟨s⟩ acts on Ω by the given action
        • ω ∈ Ω
    Out: The orbit of ω under s, i.e. ωᴳ
    =#
    ωᴳ = [x]
    γ = action(ω, s)
    while γ != x
        push!(ωᴳ, γ)
        γ = action(γ, s)
    end
    return ωᴳ
end

@inline function makeSymmetric!(S::AbstractVector{<:GroupElement})
    S = S ∪ [inv(s) for s in S]
end

function orbit(S::AbstractVector{<:GroupElement}, Ω, action=^, makeSymmetric=true)
    #=
    In:  • G = ⟨S⟩ acts on Ω by the given action
         • A set Ω
    Out: • The orbit of Ω under G, i.e. Ωᴳ
         • The transversal  # todo: what kind of transversal?  # todo: this is not correct
    Remark: If the group is infinite, the resulting orbit is only correct if S is symmetric
    =#
    @assert !isempty(S)
    if makeSymmetric
        makeSymmetric!(S)
    end
    Ωᴳ = Ω
    Ωᴳ_check = Set(Ωᴳ)
    for δ in Ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                push!(Ωᴳ, γ)
                push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ
end

function transversal(S::AbstractVector{<:GroupElement}, Ω, action=^, makeSymmetric=true)
    #=
    In:  • G = ⟨S⟩ acts on Ω by the given action
         • A set Ω
    Out: • The orbit of Ω under G, i.e. Ωᴳ
         • The transversal  # todo: what kind of transversal?
    Remark: If the group is infinite, the resulting orbit is only correct if S is symmetric
    =#
    @assert !isempty(S)
    if makeSymmetric
        makeSymmetric!(S)
    end
    Ωᴳ = Ω
    Ωᴳ_check = Set(Ωᴳ)
    T = Dict(x => one(first(S)))
    for δ in Ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                if haskey(T, δ)
					T[γ] = T[δ] * s
				end
                push!(Ωᴳ, γ)
                push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ, T
end

function transversal_factored(S::AbstractVector{<:GroupElement}, Ω, action=^, makeSymmetric=true)
    #=
    In:  • G = ⟨S⟩ acts on Ω by the given action
         • A set Ω
    Out: • The orbit of Ω under G, i.e. Ωᴳ
         • The transversal  # todo: what kind of transversal?
    Remark: If the group is infinite, the resulting orbit is only correct if S is symmetric
    =#
    @assert !isempty(S)
    if makeSymmetric
        makeSymmetric!(S)
    end
    Ωᴳ = Ω
    Ωᴳ_check = Set(Ωᴳ)
    T = Dict(x => one(first(S)))
    for δ in Ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                if haskey(T, δ)
					T[γ] = [T[δ]; s]
				end
                push!(Ωᴳ, γ)
                push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ, T
end

function schreier(S::AbstractVector{<:GroupElement}, Ω, action=^, makeSymmetric=true)  # todo: rename
    #=
    In:  • G = ⟨S⟩ acts on Ω by the given action
         • A set Ω
    Out: • The orbit of Ω under G, i.e. Ωᴳ
         • The transversal  # todo: what kind of transversal?
    Remark: If the group is infinite, the resulting orbit is only correct if S is symmetric
    =#
    @assert !isempty(S)
    if makeSymmetric
        makeSymmetric!(S)
    end
    Ωᴳ = Ω
    Ωᴳ_check = Set(Ωᴳ)
    Sch = Dict(x => one(first(S)))
    for δ in Ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                if haskey(T, δ)
					Sch[γ] = s
				end
                push!(Ωᴳ, γ)
                push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ, Sch
end

function orbit(S::AbstractVector{<:GroupElement}, ω, action=^)
    return orbit(S, [ω], action)
end


macro perm_str(s::String)
    permutation = one(Permutation([1]))
    subS = s
    while true
        m = match(r"(\((?<first>.{2,}?)\)+?)(?<rest>.*)", subS)
        firstCycle = m[:first]
        rest = m[:rest]
        parsed = parse.(Int, split(firstCycle, ","))
        permutationList = collect(1:maximum(parsed))
        for i in eachindex(parsed)
            permutationList[parsed[i]] = parsed[i % length(parsed) + 1]
        end
        permutation = :(Permutation($permutationList) * permutation)
        subS = subS[length(firstCycle) + 2: length(subS)]
        rest != "" || break
    end
    return permutation
end

function normalClosure(S::AbstractVector{<:GroupElement}, U::AbstractVector{<:GroupElement})
    N = copy(U)
    for n in N
        for s in S
            γ = inv(s) * n * s
            if γ ∉ N  # todo: it should be γ ∉ ⟨N⟩, this can be done with a sift
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
