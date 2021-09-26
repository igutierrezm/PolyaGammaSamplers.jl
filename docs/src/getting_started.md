# Getting Started

## Installation

We start by drawing 2 observations from a Polya-Gamma random variable with 
parameters 2 and 1.0.

The first step is to set up the environment:

```julia
julia> using Random, Distributions, PolyaGammaSamplers;

julia> Random.seed!(123); # Setting the seed
```

Then, we create a sampler `s`. As the first parameter is a small 
integer, we employ the PSW sampler ([1]).

```julia
julia> s = PolyaGammaPSWSampler(2, 1.0);
```

Finally, we obtain samples using `rand`.

```julia
julia> x = rand(s, 2)
5-element Vector{Float64}:
 0.19409207900581654
 0.21775655550527578
```

## References

* [1] <https://doi.org/10.1080/01621459.2013.829001>