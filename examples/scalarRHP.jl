using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots


### Scalar case
L = 1.
g = x -> .2exp(-x^2)
G = Fun(zai(g),OscLaurent(-1.,L),100) + Fun(zai(g),OscLaurent(1.,L),100)
h = x -> .5/(1 + x^2)
H = Fun(zai(h),OscLaurent(-1.,L),100)
𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)

op = x -> condense(x - G*(𝓒*x))
out = GMRES(op,H,⋅,1e-8,100, condense)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2]) ])

1.0 + (𝓢*u)(1.)
1.0 + (𝓒*u)(1.)

jump = x -> (𝓢*u)(x) - (𝓒*u)(x)*G(x) - H(x)
jump(.11)


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
