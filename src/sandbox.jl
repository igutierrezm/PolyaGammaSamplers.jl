# In the notation of [1]:
# Return K(t), given z.
function K_fun(t::Real, z::Real)
    u = t - z^2 / 2
    tol = sqrt(eps(typeof(u)))
    out = u > tol ? log(cos(sqrt(2 * u))) : log(cosh(sqrt(-2 * u)))
    out += log(cosh(z))
end

K_fun(1.0, 2.0)

# In the notation of [1]:
# Return δ(x), given x_c.
function delta_fun(x::Real, xc::Real)
    if (x <= xc)
        return (1 / xc - 1 / x) / 2
    else
        return log(x) - log(xc)
    end
end

delta_fun(1.0, 2.0)

# In the notation of [1]:
# Return δ'(x), given xc.
function delta_prime1_fun(x::Real, xc::Real)
    if (x <= xc)
        return 0.5 / x^2
    else
        return 1.0 / x
    end
end

delta_prime1_fun(1.0, 1.0)

# In the notation of [1]:
# Return φ(x), given z and u*.
function phi_fun(x::Real, z::Real, ustar::Real)
    tstar = ustar + z^2 / 2
    K_fun(tstar, z) - tstar * x
end

phi_fun(0.5, 1.0, ustar_fun(0.5))

# In the notation of [1]:
# Return φ'(x), given z and u*.
function phi_prime1_fun(x::Real, z::Real, ustar::Real)
    tstar = ustar + z^2 / 2
    -tstar
end

phi_prime1_fun(0.5, 1.0, ustar_fun(0.5))

function rel_tanh(x::Real)
    tol = typeof(x) |> eps |> sqrt
    abs(x) > tol ? tan(x) / x : 1.0 - x^2 / 3 + 2 * x^4 / 15
end

function spn_fun(x, z, ustar, n)
    Kp2star = K_prime2_star_fun(x, ustar)
    sqrt(n / 2 / pi) * exp(n * phi_fun(x, z, ustar)) / sqrt(Kp2star)
end
