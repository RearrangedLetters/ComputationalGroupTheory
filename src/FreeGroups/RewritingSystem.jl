include("Ordering.jl")

const Rule{W} = Pair{W, W} where {W <: AbstractWord}

"""
    rewrite!(v::AbstractWord, w::AbstractWord, rws::RewritingSystem)
Rewrite word `w` storing the result in `v` by left using rewriting rules of
rewriting system `rws`. See [Sims, p.66]
"""
function rewrite!(v::AbstractWord, w::AbstractWord, rule::Rule)
    left, right = rule
    while !isone(w)
        push!(v, popfirst!(w))
        if issuffix(left, v)
            prepend!(w, right)
            resize!(v, length(v) - length(left))
            break
        end
    end
    return v
end

struct RewritingSystem{O<:WordOrdering, W<:AbstractWord}
    ordering::O
    rules::AbstractVector{Rule{W}}

    function RewritingSystem(ordering::O, rules::AbstractVector{Rule{W}}) where {O<:WordOrdering, W <: AbstractWord}
        new{O, W}(ordering, rules)
    end

    function RewritingSystem(ordering::O, rules::Rule{W}...) where {O<:WordOrdering, W <: AbstractWord}
        new{O, W}(ordering, rules)
    end
end

function Base.iterate(rewritingSystem::RewritingSystem, i=1)
    return if i > length(rewritingSystem.rules) nothing else (rewritingSystem.rules[i], i + 1) end
end

function Base.empty(rws::RewritingSystem)
    return RewritingSystem(empty(rws.rwrules), ordering(rws))
end

ordering(R::RewritingSystem) = R.ordering
alphabet(lenlex::LenLex) = lenlex.alphabet
alphabet(R::RewritingSystem) = alphabet(ordering(R))
rules(R::RewritingSystem) = R.rules

function Base.show(io::IO, R::RewritingSystem)
    A = alphabet(R)
    println(io, "Rewriting system ordered by ", ordering(R), ":")
    l = ceil(Int, log10(length(rules(R))))
    ll = mapreduce(length ∘ first, max, rules(R))
    for (i, r) in enumerate(rules(R))
        println(
            io,
            " ",
            lpad(i, l),
            ". ",
            string_repr(r, A; lhspad = 2ll - 1),
        )
    end
end

function rewrite!(v::AbstractWord, w::AbstractWord, R::RewritingSystem)
    resize!(v, 0)
    while !isone(w)
        push!(v, popfirst!(w))
        for (lhs, rhs) in rules(R)
            if issuffix(lhs, v)
                prepend!(w, rhs)
                resize!(v, length(v) - length(lhs))
                break
            end
        end
    end
    return v
end

"""
Exercise: Confluent rewriting system for ℤ².
"""

function longest_common_prefix(A::AbstractWord, B::AbstractWord)
    (isone(A) || isone(B)) && return one(A)
    for i ∈ 1:min(length(A), length(B))
        if A[i] ≠ B[i]
            return A[1:i - 1]
        end
    end
    return length(A) < length(B) ? copy(A) : copy(B)
end

function isconfluent(R::RewritingSystem)
    for (P₁, Q₁) ∈ R
        for S ∈ (P₁[i:end] for i ∈ 1:length(P₁))
            for (P₂, Q₂) ∈ R
                W = longest_common_prefix(S, P₂)
                isone(W) && continue
                if length(W) == length(S)
                    A = P₁[1:length(P₁) - length(S)]
                    B = P₂[length(S) + 1:length(P₂)]
                    U = rewrite(Q₁ * B, R)
                    V = rewrite(A * Q₂, R)
                    U ≠ V && return false, (A * S * B, U, V)
                elseif length(W) == length(P₂)
                    A = P₁[1:length(P₁) - length(S)]
                    B = P₂[length(A) + length(W) + 1:length(P₁)]
                    U = rewrite(Q₁, R)
                    V = rewrite(A * Q₂, R)
                    U ≠ V && return false, (A * S * B, U, V)
                end
            end
        end
    end
    return true
end

# Example from the lecture
begin
    X = Alphabet(:x, :y, :z)
    O = LenLex(X, [:x, :y, :z])
    x, y, z = Word([X[1]]), Word([X[2]]), Word([X[3]])
    ε = one(x)
    R = RewritingSystem(O, [x * y * z => ε, y * z * x => ε, z * x * y => ε])
    isconfluent(R)
end

begin
    r = x*X => one(x)
    r isa Rule
end