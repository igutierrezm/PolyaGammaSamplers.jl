# In the notation of [1]:
# Return u*(x) := t(x) - z^2 / 2, given x.
function ustar_fun(x::Real)
    ustar0 = init_ustar(x)
    gh(u) = target_gh(u, x)
    newtons_method(gh, ustar0)
end

# In the notation of [1]:
# Return the first two derivatives of the target K(s) - s * x (wrt u).
function target_gh(u::Real, x::Real)
    gu = K_prime1_fun(u) - x
    hu = K_prime2_fun(u, gu + x)
    (gu, hu)
end

# In the notation of [1]:
# Return the 1st derivative of K(t) (wrt u)
# Note: The approximation around 0 comes from the BayesLogit (R) package
function K_prime1_fun(u::Real)
    s = 2 * u
    a = sqrt(abs(s))
    tol = sqrt(eps(typeof(s)))
    if s > tol
        return tan(a) / a
    elseif s < -tol
        return tanh(a) / a
    else
        # Use a 3rd order Taylor expansion around 0
        return 1.0 + s / 3 + 2 * s^2 / 15 + 17 * s^3 / 315
    end
end

# In the notation of [1]:
# Return the 2nd derivative of K(t) (wrt u)
# Note: The Taylor expansion around 0 comes from the BayesLogit (R) package
function K_prime2_fun(u, K_prime1)
    s = 2 * u
    if abs(s) > sqrt(eps(typeof(s)))
        return K_prime1^2 + (1 - K_prime1) / s
    else
        # Use a 1st order Taylor expansion for the second term
        return K_prime1^2 - 1 / 3 - 2 * s / 15
    end
end

# In the notation of [1]:
# Return a starting guess for u*(x) = t(x) - z^2 / 2
# Note: The guesses were taken from the (Python) package polyagamma,
# see https://github.com/zoj613/polyagamma/blob/main/src/pgm_saddle.c#L119
function init_ustar(x)
    x <= 0.25 ? -9.0 :
    x <= 0.5 ? -1.78 :
    x <= 1.0 ? -0.15 :
    x <= 1.5 ? +0.35 :
    x <= 2.5 ? +0.72 :
    x <= 4.0 ? +0.95 :
    1.15
end

# Find a root of f(x) using Newton's method, given an initial value x
# and a function `fg()` such that fg(x) returns f(x) and g(x) := f'(x).
# Note: This function was adapted from [2, p.91]
function newtons_method(fg, x; x_tol = 1e-8, iterations = 100)
    dx = Inf
    iter = 1
    fx, gx = fg(x)
    noise = gx |> typeof |> eps |> sqrt
    while abs(dx) > x_tol && iter <= iterations
        fx, gx = fg(x)
        dx = fx / (gx + noise)
        iter += 1
        x -= dx
    end
    return x
end

# tstar_fun(0.5, 1.0) # tstar(0.5, 1.0) ≈ -1.33, entonces
# tstar_fun(0.5, 1.0) - 1.0^2 / 2 # ustar(0.5) debería ser -1.83
@time ustar_fun(0.5);
