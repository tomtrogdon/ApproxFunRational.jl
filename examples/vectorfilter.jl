using  ApproxFunRational, AbstractIterativeSolvers
### Vector case
tol = 1.e-9
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
pmm = .2
L = 1.
D = 1. # overall delay

fPss = z -> 2*rs*s0s^2.0/(rs^2+z.^2)
#fPss = z -> .1*exp(-z^2/2)
fPnn = z -> pnn
fPmm = z -> pmm
Pss = Fun(zai(fPss),OscLaurent(0.0,L))
One = pad(Fun(1.,OscLaurent(0.0,L)),3)
Pnn = pnn*One
Pmm = pmm*One

H = zeros(Fun,2,2)
a = Pss + Pnn + Pmm - One
b = Pss*Fun(1.,OscLaurent(-Δs,L)) + Pnn*Fun(1.,OscLaurent(-Δn,L))
bt = Pss*Fun(1.,OscLaurent(Δs,L)) + Pnn*Fun(1.,OscLaurent(Δn,L))
H[1,1] = copy(a)
H[2,2] = copy(a)
H[1,2] = copy(b)
H[2,1] = copy(bt)
G = map(SumFun,copy(H))

# Set up right-hand side
b1 = Pss*Fun(1.,OscLaurent(- D + Ds1,L))
b2 = Pss*Fun(1.,OscLaurent(- D + Ds2,L))
h = map(SumFun,[b1,b2])

𝓒 = Cauchy(-1)
𝓢 = Cauchy(1)
simp(f) = combine!(chop!(f,tol))
inner(a,b) = ⋅(simp(a),simp(b))

function op(x::Array{T,1}) where T<:SumFun
    println("Apply Cauchy:")
    @time y = 𝓒*x
    println("Apply G:")
    @time y = G*y
    println("Subtract:")
    @time y = (x - y)
    return y
end

sop = x -> ( x - G*(𝓒*x))
out = GMRES(op,h,⋅,10*tol,30,simp)






u = sum([out[2][i]*out[1][i] for i=1:length(out[2])])

si_op(out[1][end])



𝓕 = FourierTransform(-1.0)
U = 𝓕*u

x = 0:.01:10
y1 = real(map(U[1],x))
y2 = real(map(U[2],x))

for i = 1:length(y1)-1
    if isinf(y1[i])
        y1[i] = y1[i+1]
    end
    if isinf(y2[i])
        y2[i] = y2[i+1]
    end
end

plot(x,y1)
plot!(x,y2)
