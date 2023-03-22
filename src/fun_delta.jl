# Return δ(x; xc).
function fun_δ(x::Real, xc::Real)
    x <= xc ? (1 / xc - 1 / x) / 2 : log(x) - log(xc)
end

# Return δ'(x, xc) (wrt x).
function fun_gδ(x::Real, xc::Real)
    x <= xc ? 0.5 / x^2 : 1.0 / x
end
