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
        Base.push!(ωᴳ, γ)
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
                Base.push!(Ωᴳ, γ)
                Base.push!(Ωᴳ_check, γ)
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
    for i in 1:degree(σ)
        if !visited[i]
            Δ = orbit(σ, i, ^)
            visited[Δ] .= true
            Base.push!(cycles, Δ)
        end
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

function transversal(S::AbstractVector{<:GroupElement}, x, action=^, makeSymmetric=false)
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
    Ωᴳ = [x]
    Ωᴳ_check = Set(Ωᴳ)
    T = Dict(x => one(first(S)))
    for δ in Ωᴳ
        for s in S
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                if haskey(T, δ)
					T[γ] = T[δ] * s
				end
                Base.push!(Ωᴳ, γ)
                Base.push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ, T
end

function transversalFactored(S::AbstractVector{<:GroupElement}, x, action=^, makeSymmetric=false)
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
    Ωᴳ = [x]
    Ωᴳ_check = Set(Ωᴳ)
    T = Dict(x => [one(first(S))])
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

function transversalSchreier(S::AbstractVector{<:GroupElement}, x, action=^, makeSymmetric=true)  # todo: implement
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
    Ωᴳ = [x]
    Ωᴳ_check = Set(Ωᴳ)
    schreierVector = Dict(x => 1)
    for δ in Ωᴳ
        for (i, s) in pairs(S)
            γ = action(δ, s)
            if γ ∉ Ωᴳ_check
                if haskey(schreierVector, δ)
					schreierVector[γ] = i
				end
                push!(Ωᴳ, γ)
                push!(Ωᴳ_check, γ)
            end
        end
    end
    return Ωᴳ, schreierVector
end

"""
	representative(y, S, Δ, Sch[, action=^])
Compute a representative `g` of left-coset `Stab_G(x)g` corresponding to point `y ∈ Δ` in the orbit of `x`.

## Input
* `y` - a point in `Δ`,
* `S` - a set of generators for `G = ⟨S⟩`,
* `Δ` - the orbit of `x` under the action of `G`,
* `Sch` - a Schreier tree for `Δ` and `S`.
## Output
* `r ∈ G` such that `xᵍ = y`.
"""
function representative(y, S::AbstractVector{<:GroupElement}, Δ, Sch, action=^)
    γ = y
	r = one(first(S))
	while γ ≠ first(Δ)
		s = S[Sch[γ]]
		r = s * r
		γ = action(γ, inv(s))
	end
    return r
end