"""
    PolyaGammaPSWSampler(b::Int, z::Real)

PSW sampler from [1] for a Polya-Gamma distribution with parameters `b` and `z`,
that is to say, a random variable with Laplace transform

```math
\\mathcal{L}(t) = \\cosh^b(z) \\cosh^{-b}(\\sqrt{2t + z^2})
```

References:

* [1] <https://doi.org/10.1080/01621459.2013.829001>
"""
struct PolyaGammaPSWSampler{T <: Real} <: Sampleable{Univariate, Continuous}
    b::Int
    z::T
end

function Base.rand(rng::AbstractRNG, s::PolyaGammaPSWSampler)
    out = 0.0
    for _ in 1:s.b
        out += rpg_devroye_1(rng, s.z)
    end
    return out
end

function rpg_devroye_1(rng::AbstractRNG, z::Real)
    t = 0.64
    z = abs(z) * 0.5
    yz = π^2 / 8 + z^2 / 2
    while true
        if rand(rng) < mass_texpon(z, t)
            x = t + randexp(rng) / yz
        else
            x = rtigauss(rng, z, t)
        end
        s = acoef(0, x, t)
        y = s * rand(rng)
        n = 0
        while true
            n += 1
            if n % 2 == 1
                s -= acoef(n, x, t)
                y <= s && return x / 4
            else
                s += acoef(n, x, t)
                y > s && break
            end
        end
    end
end

function mass_texpon(z::Real, t::Real)
    yz = π^2 / 8 + z^2 / 2
    b = inv(sqrt(t)) * (t * z - 1)
    a = -inv(sqrt(t)) * (t * z + 1)
    t0 = log(yz) + yz * t
    tb = t0 - z + normlogcdf(b)
    ta = t0 + z + normlogcdf(a)
    qdivp = fourinvπ * (exp(tb) + exp(ta))
    return 1.0 / (1.0 + qdivp)
end

function rtigauss(rng::AbstractRNG, z::Real, t::Real)
    z = abs(z)
    μ = inv(z)
    x = t + 1.0
    if μ > t
        α = 0.0
        while rand(rng) > α
            e1 = randexp(rng)
            e2 = randexp(rng)
            while e1^2 > (2 * e2 / t)
                e1 = randexp(rng)
                e2 = randexp(rng)
            end
            x = t / (1 + t * e1)^2
            α = exp(-0.5 * z^2 * x)
        end
    else
        while x > t
            y = randn(rng)^2
            x = μ + μ^2 * y / 2 - μ * √(4 * μ * y + (μ * y)^2) / 2
            if rand(rng) > μ / (μ + x)
                x = μ^2 / x
            end
        end
    end
    return x
end

function acoef(n::Int, x::Real, t::Real)
    if x > t
        π * (n + 0.5) * exp(-(n + 0.5)^2 * π^2 * x / 2)
    else
        2 / x / √x / sqrthalfπ * (n + 0.5) * exp(-2 * (n + 0.5)^2 / x)
    end
end

function Distributions.mean(s::PolyaGammaPSWSampler)
    s.b * inv(2.0 * s.z) * tanh(s.z / 2.0)
end

function Distributions.var(s::PolyaGammaPSWSampler)
    s.b * inv(4 * s.z^3) * (sinh(s.z) - s.z) * (sech(s.z / 2)^2)
end
