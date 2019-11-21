using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots

### Scalar case
tol = 1e-10#1.e-4
r = 1.
s0 = 1.0/sqrt(2)
pn = 1.
d = 4.0
L = 1.

Pss = z -> pn + 2*r*s0^2.0/(r^2+z.^2)
Pst = z -> 2*r*s0^2.0/(r^2+z.^2)
G = Fun(cai(Pss,pn - 1.),OscLaurent(0.0,L),100)
b = Fun(zai(Pst),OscLaurent(-d,L),100)
𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)
One = pad(Fun(1.,OscLaurent(0.0,L)),100)

op1 = x -> 𝓢*x + G*(𝓒*x)
op2 = x -> x + (G-One)*(𝓒*x)
op3 = x -> x - G*(𝓒*x) + 𝓒*x



out = GMRES(op1,b,⋅,0.000001,30, x -> x)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2])])
𝓕 = FourierTransform(-1.0)
U = 𝓕*u

x = -10:.01:10
y = real(map(U,x))
plot(x,y)

X = 0.1
U = u
causal = Cauchy(-1)*U
anti_causal = Cauchy(1)*U
(op1(u)-b)(1.)
(op2(u)-b)(1.)
(op3(u)-b)(1.)
