using ComputationalGroupTheory
using Test
using BenchmarkTools
using Profile
using PProf

S₀ = [perm"(1, 2, 3, 4)(3, 4)"]

S₁ = [perm"(1, 2, 3)", perm"(2, 3, 4)", perm"(3, 4, 1)", perm"(4, 1, 2)"]

S₂ = [perm"(1, 3, 5, 7)(2, 4, 6, 8)", perm"(1, 3, 8)(4, 5, 7)"]

S₃ = [perm"(1, 2, 3, 4)", perm"(5, 6, 7, 8)", perm"(9, 10, 11, 12)", perm"(13, 14, 15, 16)",
      perm"(1, 5, 9, 13)", perm"(2, 6, 10, 14)", perm"(3, 7, 11, 15)", perm"(4, 8, 12, 16)"]

S₄ = [perm"(2, 9)(4, 11)(6, 13)(8, 15)",
      perm"(5, 9)(6, 10)(7, 11)(8, 12)",
      perm"(3, 5)(4, 6)(11, 13)(12, 14)",
      perm"(1, 16)(2, 8)(3, 14)(4, 6)(5, 12)(7, 10)(9, 15)(11, 13)",
      perm"(1, 3)(2, 11)(4, 9)(5, 7)(6, 15)(8, 13)(10, 12)(14, 16)"]

S₅ = [perm"(1, 2, 5, 8)(3, 14, 10, 6)(4, 7, 12, 16)(9, 21, 18, 13)(11, 15, 19, 22)(17, 24, 23, 20)",
      perm"(1, 3, 5, 10)(2, 6, 8, 14)(4, 9, 12, 18)(7, 13, 16, 21)(11, 17, 19, 23)(15, 20, 22, 24)",
      perm"(1, 4, 11)(2, 7, 15)(3, 9, 17)(5, 12, 19)(6, 13, 20)(8, 16, 22)(10, 18, 23)(14, 21, 24)",
      perm"(1, 5)(2, 8)(3, 10)(4, 12)(6, 14)(7, 16)(9, 18)(11, 19)(13, 21)(15, 22)(17, 23)(20, 24)"]

S₆ = [perm"(1, 3, 8, 6)(2, 5, 7, 4)(9, 33, 25, 17)(10, 34, 26, 18)(11, 35, 27, 19)",
      perm"(9, 11, 16, 14)(10, 13, 15, 12)(1, 17, 41, 40)(4, 20, 44, 37)(6, 22, 46, 35)",
      perm"(17, 19, 24, 22)(18, 21, 23, 20)(6, 25, 43, 16)(7, 28, 42, 13)(8, 30, 41, 11)",
      perm"(25, 27, 32, 30)(26, 29, 31, 28)(3, 38, 43, 19)(5, 36, 45, 21)(8, 33, 48, 24)",
      perm"(33, 35, 40, 38)(34, 37, 39, 36)(3, 9, 46, 32)(2, 12, 47, 29)(1, 14, 48, 27)",
      perm"(41, 43, 48, 46)(42, 45, 47, 44)(14, 22, 30, 38)(15, 23, 31, 39)(16, 24, 32, 40)"]

#= @testset "SchreierSims_trivial" begin
    e = one(Permutation([1]))
    @info schreierSims([e])
    @test order([e]) == order(e)
end

@testset "SchreierSims_simple_0" begin
    σ = perm"(1, 2)"
    @info schreierSims([σ])
    @test order([σ]) == order(σ)
end =#

#= @testset "SchreierSims_S0" begin
    𝒞 = schreierSims(S₀)
    @info "Stabilizer chain: $𝒞"
    @test order(S₀) == order(first(S₀))
    @info "|G| = $(order(S₀))"
    @info "|s| = $(order(first(S₀)))"
    @test length(𝒞) == 1
end =#

#= @testset "PointStabilizer_0" begin
    𝒞 = PointStabilizer{Permutation}()
    @info 𝒞
    @test length(𝒞) == 0
end =#

#= @testset "SchreierSims_S1" begin
    @info order(S₁)
end

@testset "SchreierSims_S2" begin
    @test order(S₂) == 24
end

@testset "SchreierSims_S4" begin
    @info order(S₄)
end

@testset "SchreierSims_S5" begin
    @info order(S₅)
end

@testset "SchreierSims_S6" begin
    @info order(S₆)
end=#

#= @testset "SchreierSims_Orders" begin
    @test order(S₀) == 3
    @test order(S₁) == 12  # todo: 12 correct?
    @test order(S₃) == factorial(16)
    @test order(S₄) == 3674160  # todo: Generators probably not correct
    @test order(S₅) == 49152  # todo: correct?
    @test order(S₆) == 43252003274489856000
end  =#

#= @testset "SchreierSims_Benchmark_1" begin
    @btime order(S₆)
end =#

@testset "SchreierSims_A4" begin
    A₄ = [perm"(1, 2, 3)", perm"(2, 3, 4)"]
    println(schreierSims(A₄))
end