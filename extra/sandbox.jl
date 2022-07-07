using Revise

using BenchmarkTools
using Distributions
using PolyaGammaSamplers
using Random

# Return ``a_n(x)`` for a given t, see [1], eqs. (12)-(13)
# Algorithm 1 in [1]'s supplementary material
# Note: 
# This is a literal transcription from the article's formula 
# except for the letter case
function foo(rng, z)
    z = abs(z) / 2
    μ = 1 / z
    t = 0.64
    k = π^2 / 8 + z^2 / 2
    p = (π / 2 / k) * exp(- k * t) 
    q = 2 * exp( - z) * cdf(InverseGaussian(μ, 1.0), t)
    while true
        # Simulate a candidate x
        u = rand(rng)
        v = rand(rng)
        if (u < p / (p + q))
            # (Truncated Exponential)
            e = randexp(rng)
            x = t + e / k
        else
            # (Truncated Inverse Gaussian)
            x = randtigauss(rng, z, t)
        end
        # Evaluate if the candidate should be accepted
        s = a_xnt(x, 0, t)
        y = v * s
        n = 0
        while true
            n += 1
            if (n % 2)
                s += a_xnt(x, n, t)
                y > s && break
            else
                s -= a_xnt(x, n, t)
                y < s && return x / 4
            end
        end
    end
end

# Return ``a_n(x)`` for a given t, see [1], eqs. (12)-(13)
# Equations (12)-(139) in [1]
# Note: 
# This is a literal transcription from the article's formula
# except for the letter case
function a_xnt(x::Real, n::Int, t::Real)
    x ≤ t ? a_xnt_left(x, n, t) : a_xnt_right(x, n, t)
end

# Return ``a_n(x)^L`` for a given t
# Equation (12) in [1]
# Note: 
# This is a literal transcription from the article's formula
# except for the letter case
function a_xnt_left(x::Real, n::Int, t::Real)
    π * (n + 0.5) * (2 / π / x)^(3 / 2) * exp(- 2 * (n + 0.5)^2 / x)
end

# Return ``a_n(x)^R`` for a given t, see [1], eq. (13)
# Equation (13) in [1]
# Note: 
# This is a literal transcription from the article's formula
# except for the letter case
function a_xnt_right(x::Real, n::Int, t::Real)
    π * (n + 0.5) * exp(- (n + 0.5)^2 * π^2 * x / 2)
end

# Simulate from an IG(μ, 1) distribution
# Algorithms 2-3 in [1]'s supplementary material
# Note: 
# This is a literal transcription from the article's pseudo code
# except for the letter case
function randtigauss(rng::AbstractRNG, z::Real, t::Real)
    1 / z > t ? randtig_v1(rng, z, t) : randtig_v2(rng, z, t)
end

# Simulate from an IG(μ, 1) distribution, for μ := 1 / z > t;
# Algorithms 2 in [1]'s supplementary material
# Note:
# This is a literal transcription from the article's pseudo code
# except for the letter case and one little a detail: the 
# original condition  `x > R` must be replaced by `x > t`
function randtigauss_v1(rng::AbstractRNG, z::Real, t::Real)
    α = 0.0
    while rand(rng) > α
        e = randexp(rng) # In [1]: E 
        é = randexp(rng) # In [1]: E'
        while e^2 > (2 * é / t)
            e = randexp(rng)
            é = randexp(rng)
        end
        x = t / (1 + t * e)^2 
        α = exp(- z^2 * x / 2)
    end
    return x
end

# Simulate from an IG(μ, 1) distribution, for μ := 1 / z ≤ t
# Algorithms 3 in [1]'s supplementary material
# Note: This is a literal transcription from the article's pseudo code
function randtigauss_v2(rng::AbstractRNG, z::Real, t::Real)
    μ = 1 / z
    while x > t 
        y = randn(rng)^2
        x = μ + μ^2 * y / 2 - μ * √(4 * μ * y + (μ * y)^2) / 2
        if rand(rng) > μ / (μ + x)
            x = μ^2 / x
        end
    end
    return x
end


rng = MersenneTwister(1)
μ = 0.1
t = 2.0
z = 1 / μ
d = truncated(InverseGaussian(μ, 1.0), 0.0, t)

@benchmark PolyaGammaSamplers.rtigauss(rng, z, t)
@benchmark rand(rng, d)

# d = PolyaGammaSamplers.TruncatedIGSampler(1.0, 2.0)
# a = truncated(d, 0.0, 2.0)
# params(a)
# typeof(a)

# function foo(rng::AbstractRNG, d)
#     quantile(d, u * cdf(d, 5.0))    
# end

# foo(rng)






