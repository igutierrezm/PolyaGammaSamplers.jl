```@meta
CurrentModule = PolyaGammaSamplers
```

# PolyaGammaSamplers.jl

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
provides a method for sampling from the Polya-Gamma distribution ([1]), using 
the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) 
interface. The Polya-Gamma distribution with parameters b > 0 and z ≥ 0 
has Laplace transform

```math
\begin{aligned}
\mathcal{L}(t) = \cosh^b(z) \cosh^{-b}(\sqrt{2t + z^2}).
\end{aligned}
```

This distribution is useful for learning Bayesian Logit, 
Negative Binomial and non-parametric models (\[1, 2, 3\]).

Currently, the only available sampler is the one proposed in \[1\] (hereafter,
the PSW sampler).

## References

1. <https://doi.org/10.1080/01621459.2013.829001>.
2. <https://arxiv.org/pdf/1206.6456>.
3. <https://doi.org/10.1016/j.jspi.2020.05.009>.
4. <https://arxiv.org/abs/1405.0506>.