# Approximate u*(x) := t(x) - z^2 / 2 given x
function ustar_fun(x::Real)
    ustar0 = init_ustar(x)
    f(u) = K_prime1_fun(u) - x
    g(u) = K_prime2_fun(u)
    newtons_method(f, g, ustar0)
end

# Compute a numerically stable approximation of K'(t), given z
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

# Compute a numerically stable approximation of K''(t), given z
# Note: The Taylor expansion around 0 comes from the BayesLogit (R) package
function K_prime2_fun(u)
    K_prime1 = K_prime1_fun(u)
    s = 2 * u
    if abs(s) > sqrt(eps(typeof(s)))
        return K_prime1^2 + (1 - K_prime1) / s
    else
        # Use a 1st order Taylor expansion for the second term
        return K_prime1^2 - 1 / 3 - 2 * s / 15
    end
end

# Find a root of f(x) using Newton's method (adapted from [2, p. 91])
function newtons_method(f, g, x; x_tol = 1e-8, iterations = 1_000)
    dx = Inf
    iter = 1
    noise = g(x) |> typeof |> eps |> sqrt
    while abs(dx) > x_tol && iter <= iterations
        dx = f(x) / (g(x) + noise)
        iter += 1
        x -= dx
    end
    return x
end

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

# tstar_fun(0.5, 1.0) # tstar(0.5, 1.0) ≈ -1.33, entonces
# tstar_fun(0.5, 1.0) - 1.0^2 / 2 # ustar(0.5) debería ser -1.83
ustar_fun(0.5)

#region TODO
# - Improve get_initial_tstar().
#endregion

#region References
# [1]
# https://arxiv.org/abs/1405.0506
#endregion

#region trash
# # Compute t*(x) := \text{armgin}_s Kz(s) - s * x using bisection
# # Note: this is called t(x) in [1]
# function t_star_fun2(x::Real, z::Real)
#     f(t) = K_prime1_fun(t, z) - x
#     a = -16.0
#     b = +16.0
#     bisection(f, a, b)
# end

# # Find a root of f(x) in [a, b] (taken from [2, p.50])
# function bisection(f, a, b; x_tol = 1e-8, iterations = 1_000)
#     ya = f(a)
#     yb = f(b)
#     ya == 0 && (b = a)
#     yb == 0 && (a = b)
#     iter = 1
#     while b - a > x_tol || iter <= iterations
#         x = (a + b) / 2
#         y = f(x)
#         if y == 0
#             a, b = x, x
#         elseif sign(y) == sign(ya)
#             a = x
#         else
#             b = x
#         end
#         iter += 1
#     end
#     println(iter)
#     return (a + b) / 2
# end

# # Expand interval = [a, b] until f(a) * f(b) < 0 (taken from [2, p.50])
# function bracket_sign_change!(f, interval; k = 2)
#     a, b = interval
#     center = (b + a) / 2
#     half_width = (b - a) / 2
#     while f(a) * f(b) > 0
#         half_width *= k
#         a = center - half_width
#         b = center + half_width
#     end
#     interval[1] = a
#     interval[2] = b
#     return interval
# end
#endregion