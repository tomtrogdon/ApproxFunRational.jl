using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots


### Scalar case
L = 1.
g = x -> .4exp(-x^2)
G = Fun(zai(g),OscLaurent(1.,L),100) + Fun(zai(g),OscLaurent(0.,L),100)
h = x -> .5/(1 + x^2)
H = Fun(zai(h),OscLaurent(-1.,L),100)
𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)

op = x -> x - G*(𝓒*x)
out = GMRES(op,G,⋅,1e-10,100, x -> x)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2]) ])

l*u

l = Fun(OscConstantSpace(0.0,L),[1.])

(op(H)).space


u.space

l + u

jump = x -> (𝓢*u)(x) + 1 - ((𝓒*u)(x)+1)*(1 + G(x))
jump(.11)

jump.space

conj(jump).space

jump

𝓕 = FourierTransform(1.0)
U = 𝓕*u

x = 0:.01:10
y = real(map(U,x))
plot(x,y)




res = op(u)-b
sqrt(res⋅res)
