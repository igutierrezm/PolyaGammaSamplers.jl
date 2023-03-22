# Return φ(x, z).
function fun_φ(ustar::Real, z::Real, x::Real)
    tstar = ustar + z^2 / 2
    K_fun(tstar, z) - tstar * x
end

# Return φ'(x, z) (wrt x).
function fun_gφ(ustar::Real, z::Real)
    tstar = ustar + z^2 / 2
    -tstar
end
