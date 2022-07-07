using PolyaGammaSamplers
using Distributions
using Random
using Statistics

@testset begin
    s1 = PolyaGammaPSWSampler(2, 3.0)
    s2 = PolyaGammaPSWSampler(1, 1.0)
    @test mean(s1) ≈ 0.3017160845482888
    @test mean(s2) ≈ 0.23105857863000487
    @test var(s2) ≈ 0.03444664538852302
end

# @testset begin
#     t = 0.64
#     tol = 1e-3
#     rng = MersenneTwister(1)
#     @test abs(PolyaGammaSamplers.mass_texpon(0.0, t) - 0.5776972) < tol
#     @test abs(PolyaGammaSamplers.mass_texpon(1.0, t) - 0.4605903) < tol
#     @test abs(PolyaGammaSamplers.mass_texpon(2.0, t) - 0.2305365) < tol
#     x = [PolyaGammaSamplers.rtigauss(rng, 1.0, t) for _ in 1:10000]
#     @test abs(mean(x) - .372498) < .005
# end

@testset begin
    tol = 5
    n = 10^6
    
    s = PolyaGammaPSWSampler(1, 1.0)
    @test abs(mean(s) - mean(rand(s, n))) < tol * √(var(s) / n)
    @test abs(var(s) - var(rand(s, n))) < tol * √(2 * var(s)^2 / n)

    s = PolyaGammaPSWSampler(1, 5.0)
    @test abs(mean(s) - mean(rand(s, n))) < tol * √(var(s) / n)
    @test abs(var(s) - var(rand(s, n))) < tol * √(2 * var(s)^2 / n)
end
