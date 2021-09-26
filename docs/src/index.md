```@meta
CurrentModule = PolyaGammaSamplers
```

# PolyaGammaSamplers.jl

[PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) 
provides a method for sampling from the *Polya-Gamma distribution* ([1]), using 
the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) 
interface. The Polya-Gamma distribution with parameters _b_ and _z_ 
has Laplace transform

```math
\begin{aligned}
\mathcal{L}(t) = \cosh^b(z) \cosh^{-b}(\sqrt{2t + z^2}).
\end{aligned}
```

Although somewhat contrived, this distribution is extremely useful for learning 
Bayesian Logit models ([1]), Bayesian Negative Binomial models ([2]), and their 
generalizations (see, e.g. [3]).

Currently, the only available sampler is the one proposed in [1] (hereafter,
the PSW sampler), but I hope to add the sampler proposed in [4] during the
next months.

## References

1. Polson, N., Scott, J. & Windle, J. (2013) Bayesian inference for logistic 
    models using Pólya–Gamma latent variables, *Journal of the American 
    Statistical Association*, 108:504, 1339-1349,
    <https://doi.org/10.1080/01621459.2013.829001>.
2. Zhou, M., Li, L., Dunson, D., & Carin, L. (2012). Lognormal and Gamma 
    mixed negative Binomial regression. *Proceedings of the International 
    Conference on Machine Learning*. International Conference on Machine 
    Learning, 2012, 1343–1350.
3. Rigon, T., & Durante, D. (2021). Tractable Bayesian density regression 
    via logit stick-breaking priors. *Journal of Statistical Planning and 
    Inference*, 211, 131–142. <https://doi.org/10.1016/j.jspi.2020.05.009>.
4. Windle, J., Polson, N. G., and Scott, J. G. (2014). Sampling Polya-Gamma 
    random variates: Alternate and approximate techniques. arXiv e-prints, 
    page arXiv:1405.0506.