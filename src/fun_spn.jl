# Return sp_n(x)
function spn_fun(x::Real, z::Real, n::Real, ustar::Real)
    gQ = fun_gQ(ustar)
    hQ = fun_hQ(ustar, gQ)
    sqrt(n / 2 / pi) * exp(n * fun_phi(x, z, ustar)) / sqrt(hQ)
end
