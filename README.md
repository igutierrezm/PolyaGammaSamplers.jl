# PolyaGammaSamplers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igutierrezm.github.io/PolyaGammaSamplers.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igutierrezm.github.io/PolyaGammaSamplers.jl/dev)
[![Build Status](https://github.com/igutierrezm/PolyaGammaSamplers.jl/workflows/CI/badge.svg)](https://github.com/igutierrezm/PolyaGammaSamplers.jl/actions)
[![codecov](https://codecov.io/gh/igutierrezm/PolyaGammaSamplers.jl/branch/master/graph/badge.svg?token=yGyteQFqrS)](https://codecov.io/gh/igutierrezm/PolyaGammaSamplers.jl)

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
provides a method for sampling from the *Polya-Gamma distribution* [1], using 
the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) 
interface. See the documentation for details.

## Acknowledgments

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
is basically a Julia-1.0-compatible version of
[PolyaGammaDistribution.jl](https://github.com/currymj/PolyaGammaDistribution.jl). 
The only difference is that 
[PolyaGammaDistribution.jl](https://github.com/currymj/PolyaGammaDistribution.jl) 
defines a distribution, whereas 
[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
defines a sampler. This means that 
[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
concentrates on sampling (instead of defining the pdf, cdf and the quantiles 
of the Polya-Gamma distribution, for example). The user can learn more about 
the difference between samplers and distributions 
[here](https://juliastats.org/Distributions.jl/stable/types/).

## References

1. Polson, N., Scott, J. & Windle, J. (2013) Bayesian inference for logistic 
    models using Pólya–Gamma latent variables, *Journal of the American 
    Statistical Association*, 108:504, 1339-1349,
    <https://doi.org/10.1080/01621459.2013.829001>.