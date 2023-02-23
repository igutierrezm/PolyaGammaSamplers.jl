"""
    PolyaGammaPSWSampler(b::Int, z::Real)

PSW sampler ([1]) for a Polya-Gamma distribution. A Polya Gamma with parameters
`b` > 0 and `z` ≥ 0 has Laplace transform

```math
\\mathcal{L}(t) = \\cosh^b(z) \\cosh^{-b}(\\sqrt{2t + z^2})
```

References

* [1] <https://doi.org/10.1080/01621459.2013.829001>
"""
struct PolyaGammaPSWSampler{T <: Real} <: Sampleable{Univariate, Continuous}
    b::Int
    z::T
end

function Base.rand(rng::AbstractRNG, s::PolyaGammaPSWSampler)
    out = 0.0
    s_aux = JStarPSWSampler(s.z / 2)
    for _ in 1:s.b
        out += rand(rng, s_aux) / 4
    end
    return out
end

# function rpg_devroye_1(rng::AbstractRNG, z::Real)
#     t = 0.64
#     z = abs(z) * 0.5
#     yz = π^2 / 8 + z^2 / 2
#     while true
#         if rand(rng) < mass_texpon(z, t)
#             x = t + randexp(rng) / yz
#         else
#             x = rtigauss(rng, z, t)
#         end
#         s = acoef(0, x, t)
#         y = s * rand(rng)
#         n = 0
#         while true
#             n += 1
#             if n % 2 == 1
#                 s -= acoef(n, x, t)
#                 y <= s && return x / 4
#             else
#                 s += acoef(n, x, t)
#                 y > s && break
#             end
#         end
#     end
# end

# function mass_texpon(z::Real, t::Real)
#     yz = π^2 / 8 + z^2 / 2
#     b = inv(sqrt(t)) * (t * z - 1)
#     a = -inv(sqrt(t)) * (t * z + 1)
#     t0 = log(yz) + yz * t
#     tb = t0 - z + normlogcdf(b)
#     ta = t0 + z + normlogcdf(a)
#     qdivp = fourinvπ * (exp(tb) + exp(ta))
#     return 1.0 / (1.0 + qdivp)
# end

# function rtigauss(rng::AbstractRNG, z::Real, t::Real)
#     z = abs(z)
#     μ = inv(z)
#     x = 2 * t # any number greater than t is ok
#     if μ > t # Literal implementation of [1, algorithm 2]
#         a = 0.0 # In [1]: α
#         while rand(rng) > a
#             e1 = randexp(rng) # In [1]: E
#             e2 = randexp(rng) # In [1]: E'
#             while e1^2 > (2 * e2 / t)
#                 e1 = randexp(rng)
#                 e2 = randexp(rng)
#             end
#             x = t / (1 + t * e1)^2
#             a = exp(- z^2 * x / 2)
#         end
#     else # literal implementation of [1, algorithm 3]
#         # there is a typo in the supplement.
#         # The condition is `x > t`, not `x > R` (there is no such R)
#         while x > t
#             y = randn(rng)^2
#             x = μ + μ^2 * y / 2 - μ * √(4 * μ * y + (μ * y)^2) / 2
#             if rand(rng) > μ / (μ + x)
#                 x = μ^2 / x
#             end
#         end
#     end
#     return x
# end

# function acoef(n::Int, x::Real, t::Real)
#     if x > t
#         π * (n + 0.5) * exp(-(n + 0.5)^2 * π^2 * x / 2)
#     else
#         2 / x / √x / sqrthalfπ * (n + 0.5) * exp(-2 * (n + 0.5)^2 / x)
#     end
# end

function Distributions.mean(s::PolyaGammaPSWSampler)
    s.b * inv(2.0 * s.z) * tanh(s.z / 2.0)
end

function Distributions.var(s::PolyaGammaPSWSampler)
    s.b * inv(4 * s.z^3) * (sinh(s.z) - s.z) * (sech(s.z / 2)^2)
end

#==============================================================================#

# PSW sampler for the J* distribution [1]. The J* distribution with
# parameters 1 and `z` has Laplace transform 𝓛(t) = cos^{-z}(√2t)
struct JStarPSWSampler{T <: Real} <: Sampleable{Univariate, Continuous}
    z::T
end

# Simulate from a J*(1, z) distribution
# Algorithm 1 in [1]'s supplementary material
# Note:
# This is a literal transcription from the article's formula
# except for the letter case
function Base.rand(rng::AbstractRNG, s::JStarPSWSampler)
    z = s.z
    t = 0.64
    μ = 1 / z
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
            if (n % 2 == 1)
                s += a_xnt(x, n, t)
                y > s && break
            else
                s -= a_xnt(x, n, t)
                y < s && return x
            end
        end
    end
end

# Return ``a_n(x)`` for a given t, see [1], eqs. (12)-(13)
# Equations (12)-(13) in [1]
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

# Return a draw from IG(μ, 1) 1_{(0, t]}, see [2, pp. 3-4]
function randtigauss(rng::AbstractRNG, z::Real, t::Real)
    1 / z > t ? randtigauss_v1(rng, z, t) : randtigauss_v2(rng, z, t)
end

# Return a draw from IG(μ, 1) 1_{(0, t]}, for 1 / z > t, see [2, alg. 2]
function randtigauss_v1(rng::AbstractRNG, z::Real, t::Real)
    x = t + one(t)
    α = zero(t)
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

# Return a draw from IG(μ, 1) 1_{(0, t]}, for 1 / z ≤ t, see [2, alg. 3]
# Note: There is a typo in the document: `x > R` must be replaced by `x > t`.
function randtigauss_v2(rng::AbstractRNG, z::Real, t::Real)
    x = t + one(t)
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

#region References

# [1] https://doi.org/10.1080/01621459.2013.829001.
# [2] https://doi.org/10.1080/01621459.2013.829001 (supplementary material).

#endregion