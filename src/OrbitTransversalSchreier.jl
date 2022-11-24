include("AbstractPermutation.jl")
include("Permutation.jl")
include("CyclePermutation.jl")

@inline makeSymmetric!(s::GroupElement) = makeSymmetric!([s])

@inline function makeSymmetric!(S::AbstractVector{<:GroupElement})
    S = S ∪ [inv(s) for s in S]
end

function orbit(s::GroupElement, ω::Int, action=^)
    #=
    In: • G = ⟨s⟩ acts on Ω by the given action
        • ω ∈ Ω
    Out: The orbit of ω under s, i.e. ωᴳ
    =#
    ωᴳ = [ω]
    γ = action(ω, s)
    while γ != ω
        push!(ωᴳ, γ)
        γ = action(γ, s)
    end
    return ωᴳ
end

function orbit(S::AbstractVector{<:GroupElement}, Ω::AbstractVector, action=^, makeSymmetric=false)
    #=
    In:  • G = ⟨S⟩ acts on Ω by the given action
         • Ω a set
         • makeSymmetric signals whether inverses should be added to S
    Out: • The orbit of Ω under G, i.e. Ωᴳ
         • The transversal  # todo: what kind of transversal?  # todo: this is not correct
    Remarks:
        • If the group is infinite, the resulting orbit is only correct if S is symmetric
        • Even in the finite case it could be advantageous to make S symmetric
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

function orbit(S::AbstractVector{<:GroupElement}, ω, action=^, makeSymmetric=false)
    return orbit(S, [ω], action, makeSymmetric)
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