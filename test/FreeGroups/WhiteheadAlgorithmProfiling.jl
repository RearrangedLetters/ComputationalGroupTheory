using ComputationalGroupTheory
using BenchmarkTools
using ProfileView
using Random

Random.seed!(42)

const X = Basis(symmetric_alphabet"abc")
const Y = Basis(symmetric_alphabet"abcdef")

begin 
    @btime isprimitive_naive(word"", $X)
    @btime isprimitive_naive(word"", $Y)
    @btime isprimitive_naive(word"bc", $X)
    @btime isprimitive_naive(word"bc", $Y)
end

begin
    isprimitive_naive(word"fecd", Y)
    ProfileView.@profview isprimitive_naive(word"fecd", Y)
end

begin
    word = word"baBCBa"
    bench_naive = @benchmark isprimitive_naive($word, $Y)
    bench_nielsenfirst = @benchmark isprimitive_nielsenfirst($word, $Y)

    judge(median(bench_nielsenfirst), median(bench_naive))
end

begin
    word = word"baBCBa"
    bench_naive = @benchmark isprimitive_naive($word, $Y)
    bench_nielsenonly = @benchmark isprimitive_nielsenonly($word, $Y)

    judge(median(bench_nielsenonly), median(bench_naive))
end

begin
    randomstring = rand(collect('a':'f') ∪ collect('A':'F'), 2000)
    w = Word([Symbol(s) for s ∈ randomstring])

    @btime isprimitive($w, $Y)
end