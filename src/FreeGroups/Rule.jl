const Rule{W} = Pair{W, W} where W <: AbstractWord

begin
    r = x*X => one(x)
    r isa Rule
end

function rewrite(v::AbstractWord, w::AbstractWord, rule::Rule)
    left, right = rule
    while !isone(w)
        push!(v, popfirst!(w))
        if issuffix(left, v)
            prepend!(w, right)
            resize!(v, length(v) - length(left))
        end
    end
    return v
end

