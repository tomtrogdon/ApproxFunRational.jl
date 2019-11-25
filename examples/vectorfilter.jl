using ApproxFunOrthogonalPolynomials, ApproxFunRational,
 ApproxFunFourier, ApproxFunBase, ApproxFun, AbstractIterativeSolvers, Plots

### Vector case
tol = 1.e-4
Ds1 = 0.5
Ds2 = 1.
Δs = Ds2-Ds1
Dn1 = 1.
Dn2 = 0.5
Δn = Dn2 - Dn1
rn = 1.
s0n = 1/sqrt(2.)
rs = 1.
s0s = 1/sqrt(2.)
pnn = 1
pmm = 0.2
L = 1.
D = 1. # overall delay

fPss = z -> 2*rs*s0s^2.0/(rs^2+z.^2)
fPnn = z -> pnn
fPmm = z -> pmm
Pss = Fun(zai(fPss),OscLaurent(0.0,L),50)
One = pad(Fun(1.,OscLaurent(0.0,L)),50)
Pnn = pnn*One
Pmm = pmm*One

a = Pss + Pnn + Pmm - One
b = Pss*Fun(1.,OscLaurent(-Δs,L)) + Pnn*Fun(1.,OscLaurent(-Δn,L))
bt = Pss*Fun(1.,OscLaurent(Δs,L)) + Pnn*Fun(1.,OscLaurent(Δn,L))
H[1,1] = copy(a)
H[2,2] = copy(a)
H[1,2] = copy(b)
H[2,1] = copy(bt)
G = copy(Fun(H))

# Set up right-hand side
b1 = Pss*Fun(1.,OscLaurent(- D + Ds1,L))
b2 = Pss*Fun(1.,OscLaurent(- D + Ds2,L))
h = Fun([b1,b2])

𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)
simplify = x -> chop(condense(x),tol)

inner(a,b) = ⋅(a,b,simplify)

op = x -> simplify( x - G*(𝓒*x))
out = GMRES(op,h,inner,10*tol,20,simplify)
u = sum([out[2][i]*out[1][i] for i=1:length(out[2]) ])
u = simplify(u)

𝓕 = FourierTransform(1.0)
U = map( x->𝓕*x,Array(u))

x = 0:.011:10
y1 = real(map(U[1],x))
y2 = real(map(U[2],x))
plot(x,y1)
plot!(x,y2)


X = 0.1
U = u
causal = CauchyM(U)
anti_causal = CauchyP(U)
(anti_causal(X) - causal(X) - U(X))
(G(X) + [1. 0; 0 1])*causal(X) - anti_causal(X) - h(X)





## Test cauchy
ff = x -> exp(-x^2)
f = Fun([Fun(zai(ff),OscLaurent(1.4,1.)); Fun(zai(ff),OscLaurent(0.,1.))])
causal = CauchyM(h)
anti_causal = CauchyP(h)

anti_causal(X) - causal(X) - h(X)
flip(f::Fun{OscLaurent{DD,RR}}) where {DD,RR} = Fun(conj(f.space), conj(conj(f).coefficients))
flip(f::Fun{T}) where T<:ApproxFunBase.SumSpace = sum(map(x -> flip(x),components(f)))

H = zeros(Fun,2,2)
