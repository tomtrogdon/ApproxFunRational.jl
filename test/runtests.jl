using ApproxFunRational, ApproxFunFourier, ApproxFunBase
using Test

@testset "ApproxFunRational.jl: Evaluation and multiplication" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    @test G(.1) ≈ g(.1)*exp(α*im*.1)
    F = Fun(zai(f), OscLaurent(dom,β), 200)
    @test F(.1) ≈ f(.1)*exp(β*im*.1)
    H = F*G
    @test H(.1) ≈ F(.1)*G(.1)

    L = 2.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    @test G(.1) ≈ g(.1)*exp(α*im*.1)
    F = Fun(zai(f), OscLaurent(dom,β), 200)
    @test F(.1) ≈ f(.1)*exp(β*im*.1)
    H = F*G
    @test H(.1) ≈ F(.1)*G(.1)

    L = 0.5; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 400)
    @test G(.1) ≈ g(.1)*exp(α*im*.1)
    F = Fun(zai(f), OscLaurent(dom,β), 400)
    @test F(.1) ≈ f(.1)*exp(β*im*.1)
    H = F*G
    @test H(.1) ≈ F(.1)*G(.1)
end

@testset "ApproxFunRational.jl: Differentiation and integration" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    df = z-> (1-2z)*exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    F = Fun(zai(f), OscLaurent(dom,β), 200)
    H = F*G
    dF = F'
    dF_2 = Derivative()*F
    @test dF(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test dF_2(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test sum(H) ≈ sqrt(pi)
    @test sum(G) ≈ pi*sech(1.0*pi)
    @test sum(F) ≈ exp(-3/4 + im)*sqrt(pi)

    L = 2.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    df = z-> (1-2z)*exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    F = Fun(zai(f), OscLaurent(dom,β), 200)
    H = F*G
    dF = F'
    dF_2 = Derivative()*F
    @test dF(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test dF_2(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test sum(H) ≈ sqrt(pi)
    @test sum(G) ≈ pi*sech(1.0*pi)
    @test sum(F) ≈ exp(-3/4 + im)*sqrt(pi)

    L = 0.5; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    df = z-> (1-2z)*exp(-z^2+z)
    G = Fun(zai(g), OscLaurent(dom,α), 400)
    F = Fun(zai(f), OscLaurent(dom,β), 400)
    H = F*G
    dF = F'
    dF_2 = Derivative()*F
    @test dF(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test dF_2(.1) ≈ df(.1)*exp(β*im*.1) + (β*im)*f(.1)*exp(β*im*.1)
    @test sum(H) ≈ sqrt(pi)
    @test sum(G) ≈ pi*sech(1.0*pi)
    @test sum(F) ≈ exp(-3/4 + im)*sqrt(pi)
end
