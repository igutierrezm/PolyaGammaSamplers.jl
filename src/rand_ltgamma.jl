using Distributions
struct LeftTruncatedGammaSampler <: Sampleable{Univariate, Continuous}
    a::Float64
    b::Float64
    t::Float64
end

function Base.rand(rng::AbstractRNG, s::LeftTruncatedGammaSampler)
    (; a, b, t) = s
    if a ≈ 1.0
        return t + randexp(rng) / b
    elseif a < 1
        return t * rand_ltgamma_v1(rng, a, b * t)
    elseif b > a - 1
        return t * rand_ltgamma_v2(rng, a, b * t)
    else
        return t * rand_ltgamma_v3(rng, a, b * t)
    end
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a < 1, see [1, A4]
function rand_ltgamma_v1(rng::AbstractRNG, a::Real, b::Real)
    u = rand(rng)
    e = randexp(rng)
    x = 1.0 + e / b
    while log(u) > (a - 1.0) * log(x)
        u = rand(rng)
        e = randexp(rng)
        x = 1 + e / b
    end
    return x
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a > 1 and b > a - 1, see [1, A6]
function rand_ltgamma_v2(rng::AbstractRNG, a::Real, b::Real)
    a_minus_one = a - 1.0
    e_den = 1 - a_minus_one / b
    e = randexp(rng)
    é = randexp(rng)
    x = b + e / e_den
    while x / b - 1 + log(b / x) > é / a_minus_one
        e = randexp(rng)
        é = randexp(rng)
        x = b + e / e_den
    end
    return x / b
end

# Return a draw from G(a, b) 1_{[1, ∞)} for a > 1 and b > 0, see [1, A7]
function rand_ltgamma_v3(rng::AbstractRNG, a::Real, b::Real)
    c = (b - a + sqrt((b - a)^2 + 4 * b)) / 2 / b
    a_minus_one = a - 1
    one_minus_c = 1 - c
    x = b + randexp(rng) / c
    log_u = log(rand(rng))
    log_p = a_minus_one * log(x) - x * one_minus_c
    log_m = a_minus_one * log(a_minus_one / one_minus_c) - a_minus_one
    while log_u > log_p - log_m
        x = b + randexp(rng) / c
        log_u = log(rand(rng))
        log_p = a_minus_one * log(x) - x * one_minus_c
    end
    return x / b
end

using BenchmarkTools
rng = MersenneTwister(1)
d = truncated(Gamma(2.5, 1 / 5.5); lower = 1.0)
@btime rand_ltgamma_v1(rng, 2.5, 5.5)
@btime rand_ltgamma_v3(rng, 2.5, 5.5)
@btime rand(rng, d)

#region References
# [1] https://doi.org/10.1023/A:1018534102043
#endregion