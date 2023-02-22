# Intercept and slope of the line tangent to η(.) at x, given (z, xc, u*(x))
function tangent_to_eta(x, z, xc, ustar)
    eta_val = phi_fun(x, z, ustar) - delta_fun(x, xc)
    eta_slope = phi_prime1_fun(x, z, ustar) - delta_prime1_fun(x, xc)
    intercept = eta_val - eta_slope * x
    slope = eta_slope
    intercept, slope
end

# Example
x = 0.5
z = 1.0
xc = 2.0
ustar = ustar_fun(x)

function foo(x, z, xc, ustar)
    a, b = tangent_to_eta(x, z, xc, ustar)
    return nothing
end
@time foo(x, z, xc, ustar)