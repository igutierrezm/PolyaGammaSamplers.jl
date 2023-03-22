# Return a random draw from J*(n, z), see [1, p.23]
function rand_jstar_sp(rng::AbstractRNG, n::Real, z::Real)
    # Preprocessing
    xℓ, xc, xr = get_xℓcr(z) # support points
    L0ℓ, gLℓ = get_tangent_to_η(xℓ, z, xc) # left tangent line
    L0r, gLr = get_tangent_to_η(xr, z, xc) # right tangent line

    # Construct the unnormalized left mixture weight
    bℓ = L0ℓ # see [1, p. 21]
    pℓ = -2 * gLℓ # see [1, p. 21]
    αℓ = fun_hQ(ustarℓ, fun_gQ(ustarℓ)) / xc^3
    κℓ = exp(n / 2 / xc + n * bℓ - n * sqrt(pℓ)) / sqrt(αℓ)
    μ = 1 / sqrt(pℓ)
    λ = n
    uwℓ = κℓ * logcdf(InverseGaussian(μ, λ), xc)

    # Construct the unnormalized right mixture weight
    br = L0r
    pr = -gLr
    αr = fun_hQ(ustarr, fun_gQ(ustarr)) / xc^2
    κr = sqrt(n / 2 / pi / αr) * exp(n * br) * gamma(n) / (n * pr)^n
    uwr = κr * gamma(n, -n * gLr * xc)

    # Construct the final weights
    wl = uwl / (uwl + uwr)

    # Propose a first value
    while true
        if rand(rng) < wl
            x = n * rand_rtrunc_igauss(rng, μ / n, xc / n)
        else
            x = rand_ltrunc_gamma(rng, n, 1 / (n * pr))
        end
        u = rand(rng)
        u *= x < xc ? fun_kℓ(...) : fun_kr(...)
        if U > spn(...)
            return  n * X
        end
    end
end

# Return (xℓ(z), xc(z), xr(z)), as defined in [1, p. 23]
# References
# [1] https://arxiv.org/abs/1405.0506.
# [2] https://arxiv.org/abs/1405.0506.
function get_xℓcr(z::Real)
    xℓ = get_m(z) # see [1, p. 23]
    xc = 2.75 * xℓ # see [2]
    xr = 3.00 * xℓ # see [2]
    xℓ, xc, xr
end

# Return m(z), as defined in [1, p. 23]
# References
# [1] https://arxiv.org/abs/1405.0506.
function get_m(z::Real)
    tol = typeof(z) |> eps |> sqrt
    # Note: I'll return a 3rd order Taylor expansion if |x| ≈ 0.
    abs(x) > tol ? tanh(z) / z : 1.0 - z^2 / 3 + 2 * z^4 / 15
end

# Return the intercept and slope of the tangent line to η(x) at x,
# where η(.) is defined in [1], p. 20 (lemma 14).
# References
# [1] https://arxiv.org/abs/1405.0506.
function get_tangent_to_η(x::Real, z::Real, xc::Real)
    ustar = get_ustar(x)
    gL = get_gφ(x, z, ustar) - get_gδ(x, xc)
    L0 = get_φ(x, z, ustar) - get_δ(x, xc)
    L0, gL
end

# Return the unnormalized left weight in the mixture representation of k(.)
# References
# [1] https://arxiv.org/abs/1405.0506.
function get_uwℓ()
end

function rand_pg_sp(rng::AbstractRNG, n::Real, z::Real)
    # Construct the tangent lines
    xl = m_fun(z)
    xc = 2.75 * xl # taken from the (Python) package polyagamma, see [XX]
    xr = 3.00 * xl # taken from the (Python) package polyagamma, see [XX]
    ustarl = ustar_fun(xl)
    ustarr = ustar_fun(xr)
    L0l, gLl =  tangent_to_eta(x, z, xl, ustarl)
    L0r, gLr =  tangent_to_eta(x, z, xr, ustarr)
    # Construct the unnormalized left mixture weight
    bl = L0l
    pl = -2 * gLl
    alphal = K_prime2_star_fun(x, ustarc) / xc^3
    κl = exp(n / 2 / xc + n * bl - n * sqrt(pl)) / sqrt(alphal)
    μ = 1 / sqrt(pl)
    λ = n
    p = κl * logcdf(InverseGaussian(μ, λ), xc)
    # Construct the unnormalized right mixture weight
    br = L0r
    pr = -gLr
    alphar = K_prime2_star_fun(x, ustarc) / xc^2
    κr = sqrt(n / 2 / pi / alphar) * exp(n * br) * gamma(n) / (n * pr)^n
    q = κr * gamma(n, -n * gLr * xc)
    # Construct the final weights
    wl = p / (p + q)
    wr = q / (p + q)
    # Propose a first value
    X = rand(rng) < wl ? randtigauss(rng, z, xc) :
    U = rand(rng)
    while U >
end


# # In the notation of [1]:
# # Return K(t), given z.
# function K_fun(t::Real, z::Real)
#     u = t - z^2 / 2
#     tol = sqrt(eps(typeof(u)))
#     out = u > tol ? log(cos(sqrt(2 * u))) : log(cosh(sqrt(-2 * u)))
#     out += log(cosh(z))
# end

# K_fun(1.0, 2.0)

# # In the notation of [1]:
# # Return δ(x), given x_c.
# function delta_fun(x::Real, xc::Real)
#     if (x <= xc)
#         return (1 / xc - 1 / x) / 2
#     else
#         return log(x) - log(xc)
#     end
# end

# delta_fun(1.0, 2.0)

# # In the notation of [1]:
# # Return δ'(x), given xc.
# function delta_prime1_fun(x::Real, xc::Real)
#     if (x <= xc)
#         return 0.5 / x^2
#     else
#         return 1.0 / x
#     end
# end

# delta_prime1_fun(1.0, 1.0)

# # In the notation of [1]:
# # Return φ(x), given z and u*.
# function phi_fun(x::Real, z::Real, ustar::Real)
#     tstar = ustar + z^2 / 2
#     K_fun(tstar, z) - tstar * x
# end

# phi_fun(0.5, 1.0, ustar_fun(0.5))

# # In the notation of [1]:
# # Return φ'(x), given z and u*.
# function phi_prime1_fun(x::Real, z::Real, ustar::Real)
#     tstar = ustar + z^2 / 2
#     -tstar
# end

# phi_prime1_fun(0.5, 1.0, ustar_fun(0.5))

# function rel_tanh(x::Real)
#     tol = typeof(x) |> eps |> sqrt
#     abs(x) > tol ? tan(x) / x : 1.0 - x^2 / 3 + 2 * x^4 / 15
# end

# function spn_fun(x, z, ustar, n)
#     Kp2star = K_prime2_star_fun(x, ustar)
#     sqrt(n / 2 / pi) * exp(n * phi_fun(x, z, ustar)) / sqrt(Kp2star)
# end

# function m_fun(z::Real)
#     tol = typeof(z) |> eps |> sqrt
#     abs(z) > tol ? tanh(z) / z : 1.0 - z^2 / 3 + 2 * z^4 / 15
# end

# References
# [1] https://arxiv.org/abs/1405.0506.
# [1] https://arxiv.org/abs/1405.0506.
# [XX] https://github.com/zoj613/polyagamma/blob/main/src/pgm_saddle.c#L165.
