using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots

### Scalar case
tol = 1e-10#1.e-4
r = 1.
s0 = 1.0/sqrt(2)
pn = 1.
d = 4.0
L = 1.

Pss = z -> pn - 1. + 2*r*s0^2.0/(r^2+z.^2)
Pst = z -> 2*r*s0^2.0/(r^2+z.^2)
G = Fun(cai(Pss,pn - 1.),OscLaurent(0.0,L),100)
b = Fun(zai(Pst),OscLaurent(-d,L),100)
𝓒 = Cauchy(1)
𝓢 = Cauchy(-1)

op = x -> x - G*(𝓒*x)
out = GMRES(op,b,⋅,0.0001,50, x -> x)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2])])
𝓕 = FourierTransform(1.0)
U = 𝓕*u

x = -10:.01:10
y = real(map(U,x))
plot(x,y)
