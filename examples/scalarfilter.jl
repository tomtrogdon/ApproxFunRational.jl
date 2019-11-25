using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots

### Scalar case
tol = 1e-10#1.e-4
r = 1.
s0 = 1.0/sqrt(2)
pn = 1.
d = 4.0
L = 1.

Pss = z -> 2*r*s0^2.0/(r^2+z.^2)
#Pss = z -> exp(-z^2/2)
G = Fun(zai(Pss),OscLaurent(0.0,L),100) + Fun(OscLaurent(0.0,L),[pn])
H = Fun(zai(Pss),OscLaurent(-d,L),100)
𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)
One = pad(Fun(1.,OscLaurent(0.0,L)),100)


op = x -> x - (G-One)*(𝓒*x)

out = GMRES(op,H,⋅,0.00000001,30, x -> x)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2])])
jump = x -> (𝓢*u)(x) - (𝓒*u)(x)*G(x) - H(x)
jump(.11)

𝓕 = FourierTransform(-1.0)
U = -(𝓕*(𝓒*u))

x = -10:.01:10
y1 = real(map(U,x))
plot(x,y1)

hatH = FourierTransform(-1.0)*H
x = -10:.01:10
y1 = real(map(hatH,x))
y2 = imag(map(hatH,x))
plot(x,y1)
plot!(x,y2)



k = .3
CFun = Fun(x -> u(x)*exp(1im*k*x), -25..25,400)
sum(CFun)

sum(map(x -> x(k),components(U)))
U.space

U(.3)
