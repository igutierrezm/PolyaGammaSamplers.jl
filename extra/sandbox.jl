using Revise

using BenchmarkTools
using Distributions
using PolyaGammaSamplers
using Random

s = PolyaGammaPSWSampler(1, 1.0)
s = PolyaGammaSamplers.JStarPSWSampler(1.0)
rand(s)