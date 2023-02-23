function bounding_kernel_lf(x, z, xc, ustarc, n)
    intercept_lf, slope_lf = tangent_to_eta(x, z, xc, ustarc)
    alpha_lf = K_prime2_star_fun(x, ustarc) / xc^3
    log_out = (
        log(n / 2 / pi) / 2 -
        log(alpha_lf) / 2 +
        n / 2 / xc -
        n / 2 / x -
        3 * log(x) / 2 +
        n * (intercept_lf + slope_lf * x)
    )
    exp(log_out)
end

function bounding_kernel_rg(x, z, xc, ustarc, n)
    intercept_rg, slope_rg = tangent_to_eta(x, z, xc, ustarc)
    alpha_rg = K_prime2_star_fun(x, ustarc) / xc^2
    log_out = (
        log(n / 2 / pi) / 2 -
        log(alpha_rg) / 2 +
        n * log(xc) +
        (n - 1) * log(x) +
        n * (intercept_rg + slope_rg * x)
    )
    exp(log_out)
end


# In the notation of [1]:
# Return the 2nd derivative of K(t) (wrt u) at u*(x), given x (see [1, p.19])
function K_prime2_star_fun(x, ustar)
    x^2 + (1 - x) / 2 / ustar
end