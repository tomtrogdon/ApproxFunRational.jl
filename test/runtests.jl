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
end

@testset "ApproxFunRational.jl: Differentiation and integration" begin
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
end



j = -2
L = 1.
p = 1im*L
f = z -> ((p-z)/(z+p))^j - (-1)^j
df = z -> j*((p-z)/(p+z))^(j-1)*(-2p)/(p+z)^2
dom = PeriodicLine{false,Float64}(0.,L)

F = Fun(zai(f), OscLaurent(dom,0.), 10)
F.coefficients
dF = Derivative()*F
ApproxFunRational.fsum(dF.coefficients)

dF(.1)
df(.1)








B[:,3]
ApproxFunRational.mob_derivative_getindex(dom,1,4,2)




dF_true = Fun(zai(df), OscLaurent(dom,0.), 8)
dF_true.coefficients




dF_true(.1)
dF(.1)
