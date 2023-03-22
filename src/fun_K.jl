# Return K(t, z), as defined in [1]
function fun_K(t::Real, z::Real)
    u = t - z^2 / 2
    out = logcosh(z)
    out -= u > 0 ? log(cos(sqrt(2 * u))) : logcosh(sqrt(-2 * u))
end

# References
# [1] https://arxiv.org/abs/1405.0506.
