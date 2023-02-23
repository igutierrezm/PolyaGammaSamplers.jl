# Return a draw from a G(a, b) distribution left truncated at t, see [1].
function rand_ltgamma(rng, a, b, t)
    if a <= 1
        return t * rand_ltgamma_v1(rng, a, b * t)
    else if b > a - 1
        return t * rand_ltgamma_v2(rng, a, b * t)
    else
        return t * rand_ltgamma_v3(rng, a, b * t)
    end
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a ≤ 1, see [1, alg.4]
function rand_ltgamma_v1(rng::AbstractRNG, a::Real, b::Real)
    u = rand(rng)
    e = randexp(rng)
    x = 1 + e / b
    while u > x^(a - 1)
        u = rand(rng)
        e = randexp(rng)
        x = 1 + e / b
    end
    x
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a > 1 and b > a - 1, see [1, alg.6]
function rand_ltgamma_v2(rng::AbstractRNG, a::Real, b::Real)
    e = randexp(rng)
    é = randexp(rng)
    x = b + e / (1 - (a - 1) / b)
    while x / b - 1 + log(b / x) > é / (a - 1)
        e = randexp(rng)
        é = randexp(rng)
        x = b + e / (1 - (a - 1) / b)
    end
    x / b
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a > 1 and b > 0, see [1, alg.7]
function rand_ltgamma_v3(rng::AbstractRNG, a::Real, b::Real)
    c = (b - a + sqrt((b - a)^2 + 4 * b)) / 2 / b
    x = b + randexp(rng) / c
    u = rand(rng)
    p = x^(a - 1) * exp(-x * (1 - c))
    m = ((a - 1) / (1 - c))^(a - 1) * exp(1 - a)
    while u > p / m
        x = b + randexp(rng) / c
        u = rand(rng)
        p = x^(a - 1) * exp(-x * (1 - c))
        m = ((a - 1) / (1 - c))^(a - 1) * exp(1 - a)
    end
    x / b
end

#region References
# [1] https://doi.org/10.1023/A:1018534102043
#endregion