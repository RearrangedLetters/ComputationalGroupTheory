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

Base.empty(R::RewritingSystem) = RewritingSystem(ordering(R), empty(R.rules))
ordering(R::RewritingSystem) = R.ordering
alphabet(lenlex::LenLex) = lenlex.alphabet
alphabet(R::RewritingSystem) = alphabet(ordering(R))
rules(R::RewritingSystem) = R.rules

function string_repr(r::Rule, A::Alphabet; lhspad = 2length(first(r)) - 1, rhspad = 2length(last(r)) - 1)
    lhs, rhs = r
    L = rpad(string_repr(lhs, A), lhspad)
    R = lpad(string_repr(rhs, A), rhspad)
    return "$L → $R"
end

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

function Base.push!(R::RewritingSystem, v::AbstractWord, w::AbstractWord)
    if v == w
		return R
	end
    a = rewrite(v, R)
    b = rewrite(w, R)
    if a ≠ b
        r = b > a ? b => a : a => b
        push!(R.rules, r)
		# A = alphabet(rws)
        # @info "adding a new rule: $(string_repr(r, A))"
    # else
		# A = alphabet(rws)
		# p_str = string_repr(p, A)
        # q_str = string_repr(q, A)
        # a_str = string_repr(a, A)
        # @info "rewrites of $p_str and $q_str agree: $a_str"
    end
    return R
end

function resolve_overlaps!(R::RewritingSystem{W}, r₁::Rule, r₂::Rule) where {W}
    p₁, q₁ = r₁
    p₂, q₂ = r₂
    for s in suffixes(p₁)
        if isprefix(s, p₂)
            a = p₁[begin:end-length(s)]
            b = p₂[length(s)+1:end]
            # word a*s*b rewrites in two possible ways:
            # q₁*b and a*q₂
            # we need to resolve this local failure to confluence:
            push!(R, q₁ * b, a * q₂) # the correct rule is found in push!
		elseif isprefix(p₂, s) # i.e. p₂ is a subword in p₁
		# because rws may not be reduced
            a = p₁[begin:end-length(s)]
            b = p₁[length(a)+length(p₂)+1:end]
            # word p₁ = a*p₂*b can be rewritten in two possible ways:
            # q₁ and a*q₂*b
            push!(R, q₁, a * q₂ * b)
        end
    end
    return R
end

function reduce(R::RewritingSystem)
    @warn "reduce: not implemented yet"
    return R
end

function knuthbendix_1(R::RewritingSystem, maxrules = 100)
    S = empty(R)
    for (r₁, r₂) ∈ rules(R)
        push!(S, deepcopy(r₁), deepcopy(r₂))
    end

    for (i, r₁) in enumerate(rules(S))
        for (j,r₂) in enumerate(rules(S))
			if length(S.rules) > maxrules
                @warn "Maximum number of rules has been exceeded.
                       Try running knuthbendix with larger maxrules kwarg"
				return S
            end
			# @info (i,j)
            resolve_overlaps!(S, r₁, r₂)
            r₁ == r₂ && break
            resolve_overlaps!(S, r₂, r₁)
        end
    end
    return reduce(S)
end

#= begin
    r = x*X => one(x)
    r isa Rule
end =#

# Example 1 from the lecture
begin
    X = Alphabet(:x, :y, :z)
    O = LenLex(X, [:x, :y, :z])
    x, y, z = Word([X[1]]), Word([X[2]]), Word([X[3]])
    ε = one(x)
    R = RewritingSystem(O, [x * y * z => ε, y * z * x => ε, z * x * y => ε])
    isconfluent(R)
end

# Example 2 from the lecture
begin
    X = Alphabet(:x, :y, :z)
    O = LenLex(X, [:x, :y, :z])
    x, y, z = Word([X[1]]), Word([X[2]]), Word([X[3]])
    ε = one(x)
    R = RewritingSystem(O, [x * x => ε, y * z => ε, z * y => ε])
    isconfluent(R)
end

# Example 3 from the lecture
begin
    X = Alphabet(:x, :y, :z)
    O = LenLex(X, [:x, :y, :z])
    x, y, z = Word([X[1]]), Word([X[2]]), Word([X[3]])
    ε = one(x)
    R = RewritingSystem(O, [x * y * x * y => ε, y * x * y * x => ε])
    isconfluent(R)
end

# todo: the tests below don't work yet
let X = Alphabet([:a, :b]), O = LenLex(X, [:a, :b])
    a, b = (Word([i]) for i in 1:length(X))
    ε = one(a)
    R = RewritingSystem(O, [a^2 => ε, b^3 => ε, (a * b)^3 => ε])
    reduce(knuthbendix_1(R))
end

let X = Alphabet([:a, :b]), O = LenLex(X, [:a, :b])
    a, b = (Word([i]) for i in 1:length(X))
    ε = one(a)
    R = RewritingSystem([a^2 => ε, b^2 => ε, (a * b)^2 => ε], O)
    knuthbendix_1(R)
end

let X = Alphabet([:a, :b]), O = LenLex(X, [:a, :b])
    a, b = (Word([i]) for i in 1:length(X))
    ε = one(a)
    R = RewritingSystem([a^2 => ε, b^3 => ε, (a * b)^5 => ε], O)
    RC = reduce(knuthbendix_1(R))
end

let X = Alphabet([:a, :b, :B]), O = LenLex(X, [:a, :b, :B])
    a, b, B = (Word([i]) for i in 1:length(X))
    ε = one(a)
    R = RewritingSystem(
        [
            a^2 => ε,
            b * B => ε,
            B * b => ε,
            b^3 => ε,
            (a * b)^7 => ε,
            (a * b * a * B)^4 => ε,
        ],
        O,
    )
    RC = reduce(knuthbendix_1(R))
    @assert length(rules(RC)) == 40
    RC
end

let X = Alphabet([:a, :b, :B]), O = LenLex(X, [:a, :b, :B])
    a, b, B = (Word([i]) for i in 1:length(X))
    ε = one(a)
    R = RewritingSystem(
        [
            a^2 => ε,
            b * B => ε,
            B * b => ε,
            b^3 => ε,
            (a * b)^7 => ε,
            (a * b * a * B)^1 => ε,
        ],
        O,
    )
    RC = reduce(knuthbendix_1(R, maxrules = 200))
    RC
end