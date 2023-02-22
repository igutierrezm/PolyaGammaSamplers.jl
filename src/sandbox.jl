# n = 3
# z = 10.0
# xct = m
# xlf = 0.8 * xct
# xrg = 1.2 * xct

# # Draw proposal
# blf = 0.0
# rho_lf = -2.0
# kappa_lf = exp(3 / 2 / xct + n * blf - n * sqrt(rho_lf)) / sqrt(alf)

# Compute u := t - z^2 / 2
function u_fun(t::Real, z::Real)
    t - z^2 / 2
end

# Compute Kz(t)
function K_fun(t::Real, z::Real)
    u = u_fun(t, z)
    tol = sqrt(eps(typeof(u)))
    out = u > tol ? log(cos(sqrt(2 * u))) : log(cosh(sqrt(-2 * u)))
    out += log(cosh(z))
end

K_fun(1.0, 2.0)

# Compute a numerically stable approximation of K'(t), given z
function K_prime1_fun(t::Real, z::Real)
    s = 2 * u_fun(t, z)
    sqrt_abs_s = sqrt(abs(s))
    tol = sqrt(eps(typeof(s)))
    if s > tol
        return tan(sqrt_abs_s) / sqrt_abs_s
    elseif s < -tol
        return tanh(sqrt_abs_s) / sqrt_abs_s
    else
        # Use a 3rd order Taylor expansion around 0
        return 1.0 + s / 3 + 2 * s^2 / 15 + 17 * s^3 / 315
    end
end

@time K_prime1_fun(1.0, 2.0)

# Compute a numerically stable approximation of K''(t), given z
function K_prime2_fun(t, z)
    K_prime1 = K_prime1_fun(t, z)
    u = u_fun(t, z)
    s = 2 * u
    if abs(s) > sqrt(eps(typeof(s)))
        return K_prime1^2 + (1 - K_prime1) / s
    else
        # Use a 1st order Taylor expansion for the second term
        return K_prime1^2 - 1 / 3 - 2 * s / 15
    end
end

@time K_prime2_fun(1.0, 2.0)

# Compute \delta(x), given a value of x_c
function delta_fun(x, xc)
    if (x <= xc)
        return (1 / xc - 1 / x) / 2
    else
        return log(x) - log(xc)
    end
end

# Compute \delta'(x), given a value of x_c
function delta_prime_fun(x, xc)
    if (x <= xc)
        return 0.5 / x^2
    else
        return 1.0 / x
    end
end

function phi_fun(x, z)
    t = t_fun(x, z)
    K_fun(t, z) - t * x
end

function phi_prime_fun(x, z)
    -t_fun(x, z)
end

function tangent_to_eta_intercept(x, z, xc)
    # The tangent line has the form y = y0 + m * (x - x0), where
    # - y = η(x) = φ(x) - δ(x)
    # - m = η'(x) = φ'(x) - δ'(x)
    # To find the intercept, we only need to set x0 = 0.
    phi = phi_fun(x, z)
    delta = delta_fun(x, xc)
    phi_prime = phi_prime_fun(x, z)
    delta_prime = delta_prime_fun(x, z)
    eta_val = phi - delta
    eta_slope = phi_prime - delta_prime
    intercept = eta_val - eta_slope * x
    intercept
end

function tangent_to_eta_slope(x, z, xc)
    # The tangent line has the form y = y0 + m * (x - x0), where
    # - m = η'(x) = φ'(x) - δ'(x)
    phi_prime = phi_prime_fun(x, z)
    delta_prime = delta_prime_fun(x, xc)
    eta_slope = phi_prime - delta_prime
    slope = eta_slope
    slope
end

# Compute sp_n(x)
function sp_n(x, z, n)
    # Is my `s` a `v` for Polson?
    t = t_fun(x, z)
    phi = phi_fun(x, z)
    u = t - z^2 / 2
    s = 2 * u
    if abs(u) > sqrt(eps(typeof(u)))
        K_prime2_star = x^2 + (1 - x) / s;
    else
        K_prime2_star = x^2 - 1 / 3 - 2 * s / 15
    end
    log_out = log(0.5 * n / pi) / 2 - log(K_prime2_star) / 2 + phi / n
    exp(log_out)
end

# Compute f(t0) := K'(t0) - t0, given a value of z
function f_fun(t0, z)
    K_prime1_fun(t0, z) - t0
end

# Compute f'(t0) := K''(t0), given a value of z
function f_prime_fun(t0, z)
    K_prime2_fun(t0, z)
end
