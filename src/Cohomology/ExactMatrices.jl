using LinearAlgebra

mutable struct Subspace{T<:Any}
    B::Vector{Vector{Rational{T}}}
    P::Vector{Integer}

    function Subspace(B::Vector{Vector{Rational{T}}}) where T<:Any
        P = zeros(Integer, length(B))
        for i in 1:length(B)
            P[i] = findfirst(x -> x ≠ 0, B[i])
            for j in P[i]:length(B[i])
                B[i][j] /= B[i][P[i]]
            end
        end
        new{T}(B, P)
    end
end

Base.length(S::Subspace) = length(S.B)

function depthVector(v::Vector)
    i = findfirst(x -> x ≠ 0, v)
    return i === nothing ? (0, 0) : (i, v[i])
end

function Base.in(v::Vector{Rational{T}}, S::Subspace) where T<:Any
    c = Rational{T}[]
    for i in 1:length(S)
        a = v[S.P[i]]
        push!(c, a)
        if a ≠ 0//1
            v = v - a * S.B[i]
        end
    end
    d, a = depthVector(v)
    return d == 0
end

function Base.push!(S::Subspace, v::Vector{Rational{T}}, check=true) where T<:Any
    if !check || v ∉ S
        push!(S.B, v)
        push!(S.P, findfirst(x -> x ≠ 0, v))
    end
end

asMatrix(W::Subspace) = reduce(vcat, transpose.(W.B))

"""
Returns the reduced echelonized form of A. Row operations are used.
"""
function echelonize(A::Matrix{Rational{T}}) where T<:Any
    d, c = size(A)
    B = transpose(deepcopy(A))
    r = 1
    for t in 1:d
        for s in r:c
            if B[s, t] ≠ 0
                B[s, :] = inv(B[s, t]) * B[s, :]
                for k ∈ append!(collect(1:(r - 1)), (s + 1):c)
                    B[k, :] = B[k, :] - B[k, t] * B[s, :]
                end
                if r ≠ s
                    z = B[s, :]
                    B[s, :] = B[r, :]
                    B[r, :] = z
                end
                r += 1
                break
            end
        end
    end
    return B
end

"""
Returns subspace spanning the kernel of A, i.e. xA = 0.
"""
function nullspace(A::Matrix{Rational{T}}) where T<:Any
    B = echelonize(A)
    W = Subspace(Vector{Vector{Rational{T}}}())
    c, d = size(B)
    l = zeros(Int, c)
    r = 1
    flag = false
    for t in 1:d
        for s in r:c
            if B[s, t] ≠ 0
                l[r] = t
                r += 1
                flag = true
                break
            end
        end
        if flag
            flag = false
            continue
        end
        v = zeros(Rational{T}, d)
        v[t] = 1
        for s in 1:(r - 1)
            v[l[s]] = -B[s, t]
        end
        push!(W, v)
    end
    return W.B
end

function Base.inv(A::Matrix{Rational{T}}) where T<:Any
    d = size(A, 2)
    A⁻¹ = Matrix{Rational{T}}(I, d, d)
    flag = false
    for c ∈ 1:d
        for s ∈ c:d
            if A[s, c] ≠ 0
                M = inv(A[s, c])
                A⁻¹[s, :] = M * I[s, :]
                A[s, :] = M * A[s, :]
                for t ∈ append!(collect(1:(c - 1)), (s + 1):d)
                    A⁻¹ = A⁻¹[t, :] - A[t, c] * A⁻¹[s, :]
                    A[t] = A[t] - A[t, c] * A[s, :]
                end
                if c ≠ s
                    A⁻¹[s, :], A⁻¹[c, :] = A⁻¹[c, :], A⁻¹[s, :]
                    A[s, :], A[c, :] = A[c, :], A[s, :]
                end
                flag = true
            else
                throw("Matrix $A not invertible!")
            end
            if flag
                flag = false
                break
            end
        end
    end
    return A⁻¹
end

function rowspace(A::Matrix{Rational{T}}) where T<:Any
    A = echelonize(Base.copy(A))
    l = 0
    for r ∈ 1:(size(A)[1])
        if isnothing(findfirst(x -> x ≠ 0, A[r, :])) l = r, break end
    end
    return A[1:(l - 1), :]
end
