using ComputationalGroupTheory
using Test

@testset "Count Free Group Automorphisms" begin
    #=
    First we assert that the iteration protocol actually does the desired number
    of iterations. This number is not equal to the actual number of Whitehead
    automorphisms because the Nielsen automorphisms are covered twice.
    =#
    X = symmetric_alphabet"a"
    number_of_automorphisms = length(WhiteheadAutomorphisms(X, word"a"))
    @test number_of_automorphisms == 2

    v = word"a"

    for wordlength ∈ 2:5
        word = v^wordlength
        W = WhiteheadAutomorphisms(X, word)
        @test length(W) == factorial(wordlength) * big(2)^wordlength
    end

    #=
    Now we check if we obtain the correct number of Whitehead automorphisms.
    According to [HAR] there are 
    =#
end

#= @testset "Count Nielsen Automorphisms" begin
    v = word"a"

    for wordlength ∈ 1:6
        word = v^wordlength
        Nᵢ = NielsenAutomorphisms(rewritingSystem, word)
        length_Nᵢ = 0
        for _ in Nᵢ length_Nᵢ += 1 end
        @test length_Nᵢ == 5 * wordlength^2
    end
end

@testset "Primitive elements in ℤ" begin
    X = Alphabet(:𝟙)
    setinverse!(X, :𝟙, :𝟙⁻)

    @test freeRewriteBV!(word"𝟙𝟙⁻", X) == word""

    #=
    Now F₁ ≅ ⟨𝟙⟩ ≅ ℤ, and there are exactly two primitive elements,
    namely 𝟙 ≙ 1 and -1 ≙ 𝟙⁻.
    =#
    @test isprimitive_naive(X, word"𝟙")
    @test isprimitive_naive(X, word"𝟙⁻")

    #=
    Now we assert that neither 2 ≙ 𝟙𝟙 nor -2 ≙ 𝟙⁻𝟙⁻ are primitive.
    =#
    @test !isprimitive_naive(X, word"𝟙𝟙")
    @test !isprimitive_naive(X, word"𝟙⁻𝟙⁻")
end =#