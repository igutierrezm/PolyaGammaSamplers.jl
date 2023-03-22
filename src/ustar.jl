# Return u*(x) := t(x) - z^2 / 2, with t(x) defined in [1]
function ustar(x::Real)
    ustar0 = init_ustar(x)
    gh(u) = target_gh(u, x)
    newtons_method(gh, ustar0)
end

# First 2 derivatives of Q(u) - u * x (wrt u).
function target_gh(u::Real, x::Real)
    g = gQ(u)
    h = hQ(u, g)
    g - x, h
end

# Return Q'(u).
function gQ(u::Real)
    s = 2 * u
    a = sqrt(abs(s))
    tol = sqrt(eps(typeof(a)))
    if a > tol
        return tan(a) / a
    elseif a < -tol
        return tanh(a) / a
    else
        # Use a 3rd order Taylor expansion around 0
        return 1.0 + s / 3 + 2 * s^2 / 15 + 17 * s^3 / 315
    end
end

# Return Q''(u).
function hQ(u::Real, gQ::Real)
    s = 2 * u
    if abs(s) > sqrt(eps(typeof(s)))
        return gQ^2 + (1 - gQ) / s
    else
        # Use a 1st order Taylor expansion for the second term
        return gQ^2 - 1 / 3 - 2 * s / 15
    end
end

# Return a starting guess for u*(x).
# Note: The guesses were taken from [2] (in the Python package `polyagamma`).
function init_ustar(x)
    x <= 0.25 ? -9.0 :
    x <= 0.5 ? -1.78 :
    x <= 1.0 ? -0.15 :
    x <= 1.5 ? +0.35 :
    x <= 2.5 ? +0.72 :
    x <= 4.0 ? +0.95 :
    1.15
end

# Find a root of ``f(x)`` using Newton's method, given an initial value
# ``x`` and a function `fg()` such that `fg(x)` returns ``f(x)`` and ``f'(x)``.
# Note: This function was adapted from [3, p. 91] (Algorithms for Optimization).
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

# References
# [1] https://arxiv.org/abs/1405.0506.
# [2] https://github.com/zoj613/polyagamma/blob/main/src/pgm_saddle.c#L119.
# [3] https://mitpress.mit.edu/9780262039420/algorithms-for-optimization/.