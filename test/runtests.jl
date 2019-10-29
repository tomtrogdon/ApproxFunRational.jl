using ApproxFunRational, ApproxFunFourier, ApproxFunBase
using Test

@testset "ApproxFunRational.jl: Oscillatory Cauchy integrals" begin
    L = 1.; α = -2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    GCP = CauchyP(G)
    GCM = CauchyM(G)
    @test GCP(.1) - GCM(.1) ≈ G(.1)
    @test GCP(.1) ≈ 0.027382689548799077 + 0.0017439616653601026im
    @test GCM(.1) ≈ -0.9478038907592183 + 0.1994240679809568im

    L = 2.; α = -2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    G = Fun(zai(g), OscLaurent(dom,α), 300)
    GCP = CauchyP(G)
    GCM = CauchyM(G)
    @test GCP(.1) - GCM(.1) ≈ G(.1)
    @test GCP(.1) ≈ 0.027382689548799077 + 0.0017439616653601026im
    @test GCM(.1) ≈ -0.9478038907592183 + 0.1994240679809568im

    L = 0.5; α = -2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    G = Fun(zai(g), OscLaurent(dom,α), 400)
    GCP = CauchyP(G)
    GCM = CauchyM(G)
    @test GCP(.1) - GCM(.1) ≈ G(.1)
    @test GCP(.1) ≈ 0.027382689548799077 + 0.0017439616653601026im
    @test GCM(.1) ≈ -0.9478038907592183 + 0.1994240679809568im
end

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

    f = z -> exp(-z^2+z*1im)*(1.0+1im)
    F = Fun(zai(f), OscLaurent(dom,β), 200)
    Fbar = conj(F)
    @test Fbar(.1) ≈ conj(F(.1))
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
