using ApproxFunRational, ApproxFunFourier, ApproxFunBase
using Test

@testset "ApproxFunRational.jl" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    ff = zai(f)
    gg = zai(g)
    G = Fun(gg, OscLaurent(dom,α), 600)
    @test G(.1) ≈ g(.1)*exp(α*im*.1)
    F = Fun(ff, OscLaurent(dom,β), 600)
    @test F(.1) ≈ f(.1)*exp(β*im*.1)
    H = F*G
    @test H(.1) ≈ F(.1)*G(.1)
end

L = .9; α = -2.; β = 2.;
dom = PeriodicLine{false,Float64}(0.,L)
g = z -> sech(z)
f = z -> exp(-z^2+z)
ff = zai(f)
gg = zai(g)
G = Fun(gg, OscLaurent(dom,α), 200)
F = Fun(ff, OscLaurent(dom,β), 200)
H = F*G

F(.1)*G(.1)
H(.1)
H.coefficients

S(a) = ApproxFunRational.fsum(a.coefficients)
H1 = Fun(OscLaurent(dom,0.),S(F)*G.coefficients)
H2 = Fun(OscLaurent(dom,0.),S(G)*F.coefficients)
H(.1) - H1(.1) - H2(.1)
F(.1)*G(.1)
