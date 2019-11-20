using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots

### Scalar case
L = 1.
g = x -> .2exp(-x^2)
G1 = Fun(zai(g),OscLaurent(1.,L),100) + Fun(zai(g),OscLaurent(0.,L),100)
h = x -> .2/(1 + x^2)
H2 = Fun(zai(h),OscLaurent(-1.,L),100)
𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)

G=zeros(Fun,2,2)
G[1,1] = G1
G[2,2] = G1
G[1,2] = H2
G[2,1] = H2
G = copy(Fun(G))

One = pad(Fun(1.,OscLaurent(0.0,L)),40)
b = G*Fun([One;0*One])

op = x -> condense(x - G*(𝓒*x))
out = GMRES(op,b,⋅,1e-3,100, condense)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2]) ])

id = Array([1. 0.; 0. 1.])
𝓢 = Cauchy(1)
𝓒 = Cauchy(-1)
c = id*[1.0; 0.]

id + G(.1)

jump = x -> (𝓢*u)(x) + c - (id + G(x))*((𝓒*u)(x)+c)
jump(.11)
