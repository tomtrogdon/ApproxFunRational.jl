using ApproxFunRational, ApproxFunFourier, ApproxFunBase
using Test

@testset "ApproxFunRational.jl" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    ff = zai(f)
    gg = zai(g)
    G = Fun(gg, OscLaurent(dom,α), 200)
    @test G(.1) ≈ g(.1)*exp(α*im*.1)
    F = Fun(ff, OscLaurent(dom,β), 200)
    @test F(.1) ≈ f(.1)*exp(β*im*.1)
    H = F*G
    @test H(.1) ≈ F(.1)*G(.1)
    sum(H) ≈ sqrt(pi)
    sum(G) ≈ pi*sech(1.0*pi)
    sum(F) ≈ exp(-3/4 + im)*sqrt(pi)
end
