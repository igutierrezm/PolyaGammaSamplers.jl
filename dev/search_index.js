var documenterSearchIndex = {"docs":
[{"location":"getting_started/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"getting_started/#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"We start by drawing 2 observations from a Polya-Gamma random variable with  parameters 2 and 1.0.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"The first step is to set up the environment:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> using Random, Distributions, PolyaGammaSamplers;\n\njulia> Random.seed!(123); # Setting the seed","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Then, we create a sampler s. As the first parameter is a small  integer, we employ the PSW sampler ([1]).","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> s = PolyaGammaPSWSampler(2, 1.0);","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Finally, we obtain samples using rand.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> x = rand(s, 2)\n5-element Vector{Float64}:\n 0.19409207900581654\n 0.21775655550527578","category":"page"},{"location":"getting_started/#References","page":"Getting Started","title":"References","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"[1] https://doi.org/10.1080/01621459.2013.829001","category":"page"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Install with the Julia package manager Pkg:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"# Press ']' to enter the Pkg REPL mode.\npkg> add PolyaGammaSamplers","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"julia> using Pkg; \njulia> Pkg.add(\"PolyaGammaSamplers\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = PolyaGammaSamplers","category":"page"},{"location":"#PolyaGammaSamplers.jl","page":"Home","title":"PolyaGammaSamplers.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PolyaGammaSamplers.jl  provides a method for sampling from the Polya-Gamma distribution ([1]), using  the Distributions.jl  interface. The Polya-Gamma distribution with parameters b and z  has Laplace transform","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nmathcalL(t) = cosh^b(z) cosh^-b(sqrt2t + z^2)\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"Although somewhat contrived, this distribution is extremely useful for learning  Bayesian Logit models ([1]), Bayesian Negative Binomial models ([2]), and  their generalizations (see, e.g. [3]).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, the only available sampler is the one proposed in [1] (hereafter, the PSW sampler), but I hope to add the sampler proposed in [4] soon.","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Polson, N., Scott, J. & Windle, J. (2013) Bayesian inference for logistic   models using Pólya–Gamma latent variables, Journal of the American   Statistical Association, 108:504, 1339-1349,  https://doi.org/10.1080/01621459.2013.829001.\nZhou, M., Li, L., Dunson, D., & Carin, L. (2012). Lognormal and Gamma   mixed negative Binomial regression. Proceedings of the International   Conference on Machine Learning. International Conference on Machine   Learning, 2012, 1343–1350.\nRigon, T., & Durante, D. (2021). Tractable Bayesian density regression   via logit stick-breaking priors. Journal of Statistical Planning and   Inference, 211, 131–142. https://doi.org/10.1016/j.jspi.2020.05.009.\nWindle, J., Polson, N. G., and Scott, J. G. (2014). Sampling Polya-Gamma   random variates: Alternate and approximate techniques. arXiv e-prints,   page arXiv:1405.0506.","category":"page"},{"location":"samplers/#Samplers","page":"Samplers","title":"Samplers","text":"","category":"section"},{"location":"samplers/","page":"Samplers","title":"Samplers","text":"Modules = [PolyaGammaSamplers]","category":"page"},{"location":"samplers/#PolyaGammaSamplers.PolyaGammaPSWSampler","page":"Samplers","title":"PolyaGammaSamplers.PolyaGammaPSWSampler","text":"PolyaGammaPSWSampler(b::Int, z::Real)\n\nPSW sampler ([1]) for a Polya-Gamma distribution with parameters b and z, and Laplace transform\n\nmathcalL(t) = cosh^b(z) cosh^-b(sqrt2t + z^2)\n\nReferences\n\n[1] https://doi.org/10.1080/01621459.2013.829001\n\n\n\n\n\n","category":"type"}]
}
