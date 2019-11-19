using ApproxFunRational, ApproxFunFourier, ApproxFunBase, ApproxFun
using Test

## TODO: test inner product

@testset "ApproxFunRational.jl: Fourier transform" begin
    f = x -> 1. + (x+x^2)*exp(-x^2+1im*x^2)
    α = 2.0
    L = 2.0
    F = Fun(cai(f,1.0),OscLaurent(α,L))
    𝓕 = FourierTransform(3.0)
    𝓗 = FourierTransform(-3.0)
    FF = 𝓗*(𝓕*F)
    @test F(.1) ≈ FF(.1)

    f = x -> (x+x^2)*exp(-x^2+1im*x^2)
    α = 2.0
    L = 2.0
    F = Fun(zai(f),OscLaurent(α,L))
    𝓕 = FourierTransform(1.0)
    Fhat = 𝓕*F
    @test Fhat(.1) ≈ -0.1979914917932347 + 0.2718859895675688im
end

@testset "ApproxFunRational.jl: Matrix-vector function product" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    h = z -> 1.0 + 1/(z^2+1)
    F = z -> [cai(h,1.0)(z) zai(g)(z); zai(f)(z) cai(h,1.0)(z)]
    G = z -> [zai(g)(z), zai(f)(z)]
    FF = Fun(F,OscLaurent(0.0))
    GG = Fun(G,OscLaurent(0.0))
    J = z -> [zai(g)(z)^2 + zai(f)(z)^2]
    JJ = transpose(GG)*GG
    HH = FF*GG
    @test JJ(.1) ≈ J(.1)
    @test HH(.1) ≈ F(.1)*G(.1)
    FF = Fun(F,OscLaurent(β))
    GG = Fun(G,OscLaurent(0.0)) + Fun(G,OscLaurent(α))
    @test (FF*GG)(.1) ≈ F(.1)*G(.1)*exp(1im*β*(.1)) + F(.1)*G(.1)*exp(1im*(α+β)*(.1))
end

@testset "ApproxFunRational.jl: Vector-valued functions" begin
    L = 1.; α = -2.; β = 2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    f = z -> exp(-z^2+z)
    h = z -> 1.0 + 1/(z^2+1)
    F = z-> [zai(g)(z), zai(f)(z), cai(h,1.0)(z)]
    FF = Fun(F,OscLaurent(α))
    @test F(.1)*exp(1im*α*(.1)) ≈ FF(.1)
end

@testset "ApproxFunRational.jl: Oscillatory Cauchy integrals" begin
    L = 1.; α = -2.;
    dom = PeriodicLine{false,Float64}(0.,L)
    g = z -> sech(z)
    G = Fun(zai(g), OscLaurent(dom,α), 200)
    GCP = ApproxFunRational.CauchyP(G)
    GCM = ApproxFunRational.CauchyM(G)
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
    G = Fun(zai(g), OscLaurent(dom,α))
    println(length(G.coefficients))
    @test G(.1) ≈ g(.1)*exp(α*im*.1) #test adaptivity

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
    G = Fun(zai(g), OscLaurent(dom,α), 500)
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
    G = Fun(zai(g), OscLaurent(dom,α), 400)
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
