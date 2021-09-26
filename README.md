# PolyaGammaSamplers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igutierrezm.github.io/PolyaGammaSamplers.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igutierrezm.github.io/PolyaGammaSamplers.jl/dev)
[![Build Status](https://github.com/igutierrezm/PolyaGammaSamplers.jl/workflows/CI/badge.svg)](https://github.com/igutierrezm/PolyaGammaSamplers.jl/actions)
[![codecov](https://codecov.io/gh/igutierrezm/PolyaGammaSamplers.jl/branch/master/graph/badge.svg?token=yGyteQFqrS)](https://codecov.io/gh/igutierrezm/PolyaGammaSamplers.jl)

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
provides a method for sampling from the *Polya-Gamma distribution*, as
described in [1], using the 
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl) 
interface. See the documentation for details.

# Acknowledgment

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
is basically a Julia-1.0-compatible version of
[PolyaGammaDistribution.jl](https://github.com/currymj/PolyaGammaDistribution.jl). 
The only difference is that `PolyaGammaDistribution.jl` defines a distribution,
whereas `PolyaGammaSamplers.jl` defines a sampler. This means that 
`PolyaGammaSamplers.jl` concentrates on sampling (instead of 
defining the pdf, cdf and the quantiles of the Polya-Gamma 
distribution, for example). 