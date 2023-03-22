module PolyaGammaSamplers

using Distributions
using LogExpFunctions: logcosh
using Random
using SpecialFunctions
using StatsFuns

export PolyaGammaPSWSampler

include("polyagammapswsampler.jl")

end
