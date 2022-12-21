#= 
Notes:
    • Currently the commutative ring K is not taken into account
    • Relators wᵢ ∈ R are assumed to be arrays of tuples (xᵢ, εᵢ)
      where xᵢ ∈ X and ε = ±1. For example
      S₃ ≅ ⟨ X = {f, g} | R = {[f³ = [(f, 1), (f, 1), (f, 1)], g² = [(g, 1), (g, 1)],
                            gfgf = [(g, 1), (f, 1), (g, 1), (f, 1)]]} ⟩
=#

function Z1(X, R, Λ)
    r = length(X)
    s = length(R)
    d = size(Λ[1], 2)
    Z = zeros(Rational{BigInt}, d * r, d * s)
    for k ∈ 1:s
        c = d * (k - 1)
        α = Matrix{BigInt}(I, d, d)
        wₖ = R[k]
        for i ∈ length(wₖ):-1:1
            e = wₖ[i]
            if e == -1
                α = inv(Λ[j]) * α
            end
            r = d * (j - 1)
            for i₁ ∈ 1:d
                for i₂ ∈ 1:d
                    Z[r + i₁, c + i₂] += e * α[i₁, i₂]
                end
            end
            if e == 1
                α = Λ[j] * α
            end
        end
    end
    return Z
end

function B1(X, R, Λ)
    r = length(X)
    s = length(R)
    d = size(Λ[1], 2)
    B = zeros(Rational{BigInt}, d * r, d * s)
    for j ∈ 1:r
        c = d * (j - 1)
        for i₁ ∈ 1:d
            for i₂ ∈ 1:d
                if i₁ == i₂
                    B[i₁, c + i₂] = 1 - Λ[i₁, i₂]
                else
                    B[i₁, c + i₂] = -Λ[i₁, i₂]
                end
            end
        end
    end
    return B
end