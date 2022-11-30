include("AbstractPermutation.jl")

struct Permutation <: AbstractPermutation
    images::Vector{Int}

    function Permutation(v::Vector{<:Integer}, check=true)
        if check
            @assert sort(v) == 1:length(v) "Image vector doesn't define a permutation!"
        end
        return new(v)
    end
end

function Base.:^(n::Integer, σ::Permutation)
    if n > length(σ.images)
        return convert(Int, n)
    else
        return σ.images[n]
    end
end

function Base.:^(σ::Permutation, n::Integer)
    τ = one(σ)
    k = n // 1
    while k ≥ 2
        if iseven(k)
            σ *= σ
            k //= 2
        else
            τ *= σ
            σ *= σ
            k = (k - 1) // 2
        end
    end
    return σ * τ
end

function Base.:^(σ::Permutation, τ::Permutation)
    return σ * τ
end

function Base.:(*)(σ::Permutation, τ::Permutation)
    range = max(degree(σ), degree(τ))
    result = collect(1:range)
    for i in 1:range
        result[i] = (i^σ)^τ
    end
    return Permutation(result)
end

function degree(σ::Permutation)
    return something(findlast(i -> σ.images[i] != i, 1:length(σ.images)), 1)
end

macro perm_str(s::String)  # todo: this can be optimized by calculating the images and calling the constructor only once
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
        permutation = :(Permutation($permutationList) * $permutation)
        subS = subS[length(firstCycle) + 2: length(subS)]
        rest != "" || break
    end
    return permutation
end