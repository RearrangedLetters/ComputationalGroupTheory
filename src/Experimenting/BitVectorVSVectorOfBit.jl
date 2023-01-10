using BenchmarkTools, LoopVectorization

function variant1(a, b)
    s = length(a)
    c = Vector{Bool}(undef, s)
    @inbounds @simd for i in 1:s
        c[i] = a[i] > b[i]
    end
    c
end

function variant2(a, b)
    s = length(a)
    c = BitVector(undef, s)
    @avx for i in 1:s
        c[i] = a[i] > b[i]
    end
    c
end

const a = rand(100_000)
const b = rand(100_000)

@btime variant1($a, $b)
@btime variant2($a, $b)
