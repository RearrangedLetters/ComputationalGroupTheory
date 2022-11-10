struct Permutation
    images::Vector{Int}

    function Permutation(v::Vector{<:Integer}, check=false)
        if check
            @assert sort(v) == 1:length(v) "Image vector doesn't define a permutation."
        end
        return new(v)
    end
end

(σ::Permutation)(n::Integer) = n > length(σ.images) ? convert(Int, n) : σ.images[n]

degree(σ::Permutation) =
    something(findlast(i -> σ.images[i] != i, 1:length(σ.images)), 1)

function Base.:(==)(σ::Permutation, τ::Permutation)
    deg = max(degree(σ), degree(τ))
    for i in 1:deg
        if σ(i) != τ(i)
            return false
        end
    end
    return true
end

Base.one(σ::Permutation) = Permutation([1])

function Base.inv(σ::Permutation)
    deg = degree(σ)
    τ = collect(1:deg)
    for i in deg
        τ[σ(i)] = i
    end
    return Permutation(τ)
end

function Base.:(*)(σ::Permutation, τ::Permutation)
    range = max(degree(σ), degree(τ))
    result = collect(1:range)
    for i in 1:range
        result[i] = τ(σ(i))
    end
    return Permutation(result)
end

function Base.:^(i::Integer, σ::Permutation)
    τ = one(σ)
    for _ in 1:i
        τ *= τ
    end
    return τ(i)
end

function cycleDecomposition(σ::Permutation)
    deg = degree(σ)
    visited = falses(deg)
    cycles = Vector{Vector{Int}}()
    for i in 1:deg
        if !visited[i]
            Δ = orbit_plain(i, σ)
            visited[Δ] .= true
            push!(cycles, Δ)
        end
    end
    return cycles
end

function Base.show(io::IO, σ::Permutation)
    if isone(σ)
        print(io, "()")
    end

    for cycle in cycleDecomposition(σ)
        if length(cycle) > 1
            print(io, "(")
            join(io, cycle, ",")
            print(io, ")")
        end
    end
end

function orbit_vanilla(S::AbstractVector{<:Permutation}, ω, action=^)
    #=
    G = ⟨S⟩ acts on Ω by the given action
    ω ∈ Ω
    =#
    @assert !isempty(S)
    ωᴳ = [ω]
    for δ in ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ ωᴳ
                push!(ωᴳ, γ)
            end
        end
    end
    return ωᴳ
end

# orbit_plain(x, s, action=^) = orbit_plain!([x], [s], action)
# orbit_plain(x, S::AbstractVector{<:Permutation}, action=^) = orbit_plain!([x], S, action)