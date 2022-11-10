@testset "Orbit" begin
    using ComputationalGroupTheory
    σ = Permutation([2,3,4,1])
    τ = Permutation([2,1])
    S = [σ, τ]
    Δ = orbit_vanilla(S, one(σ), *)
    @test Δ == unique(Δ)
    @test length(Δ) == 24
    @test σ * τ in Δ
end