module ApproxFunRational
using Base, ApproxFun, ApproxFunBase, ApproxFunFourier, Reexport, FFTW, LinearAlgebra, ApproxFunOrthogonalPolynomials, ApproxFunRational#, MacroTools#, Reexport, AbstractFFTs, FFTW, InfiniteArrays, FillArrays, FastTransforms, IntervalSets,
            #DomainSets, SpecialFunctions

@reexport using ApproxFunBase

#import ApproxFunBase: Fun, AbstractTransformPlan, TransformPlan, ITransformPlan, ConcreteDerivative, ConcreteMultiplication, prectype, spacescompatible

import ApproxFunBase: normalize!, flipsign, FiniteRange, Fun, MatrixFun, UnsetSpace, VFun, RowVector,
                UnivariateSpace, AmbiguousSpace, SumSpace, SubSpace, WeightSpace, NoSpace, Space,
                HeavisideSpace, PointSpace, dimension, colstop, israggedbelow, rowstop,
                IntervalOrSegment, RaggedMatrix, AlmostBandedMatrix, chop, chop!,
                AnyDomain, ZeroSpace, ArraySpace, TrivialInterlacer, BlockInterlacer,
                AbstractTransformPlan, TransformPlan, ITransformPlan,
                ConcreteConversion, ConcreteMultiplication, ConcreteDerivative, ConcreteIntegral, CalculusOperator,
                ConcreteVolterra, Volterra, VolterraWrapper,
                MultiplicationWrapper, ConversionWrapper, DerivativeWrapper, Evaluation, EvaluationWrapper,
                Conversion, defaultConversion, defaultcoefficients, default_Fun, Multiplication, Derivative, Integral, bandwidths,
                ConcreteEvaluation, ConcreteDefiniteLineIntegral, ConcreteDefiniteIntegral, ConcreteIntegral,
                DefiniteLineIntegral, DefiniteIntegral, ConcreteDefiniteIntegral, ConcreteDefiniteLineIntegral, IntegralWrapper,
                ReverseOrientation, ReverseOrientationWrapper, ReverseWrapper, Reverse, NegateEven,
                Dirichlet, ConcreteDirichlet, DirichletWrapper,
                TridiagonalOperator, SubOperator, Space, @containsconstants, spacescompatible,
                hasfasttransform, canonicalspace, domain, setdomain, prectype, domainscompatible,
                plan_transform, plan_itransform, plan_transform!, plan_itransform!, transform, itransform, hasfasttransform,
                CanonicalTransformPlan, ICanonicalTransformPlan,
                Integral, ∞,
                domainspace, rangespace, boundary,
                union_rule, conversion_rule, maxspace_rule, conversion_type, maxspace, hasconversion, points,
                rdirichlet, ldirichlet, lneumann, rneumann, ivp, bvp,
                linesum, differentiate, integrate, linebilinearform, bilinearform,
                UnsetNumber, coefficienttimes, subspace_coefficients, sumspacecoefficients, specialfunctionnormalizationpoint,
                Segment, IntervalOrSegmentDomain, PiecewiseSegment, isambiguous, Vec, eps, isperiodic,
                arclength, complexlength,
                invfromcanonicalD, fromcanonical, tocanonical, fromcanonicalD, tocanonicalD, canonicaldomain, setcanonicaldomain, mappoint,
                reverseorientation, checkpoints, evaluate, mul_coefficients, coefficients, coefficientmatrix, isconvertible,
                clenshaw, ClenshawPlan, sineshaw,
                toeplitz_getindex, toeplitz_axpy!, sym_toeplitz_axpy!, hankel_axpy!, ToeplitzOperator, SymToeplitzOperator, hankel_getindex,
                SpaceOperator, ZeroOperator, InterlaceOperator,
                interlace!, reverseeven!, negateeven!, cfstype, pad!, alternatesign!, mobius,
                extremal_args, hesseneigvals, chebyshev_clenshaw, recA, recB, recC, roots,splitatroots,
                chebmult_getindex, intpow, alternatingsum,
                domaintype, diagindshift, rangetype, weight, isapproxinteger, default_Dirichlet, scal!, dotu,
                components, promoterangespace, promotedomainspace, choosedomainspace,
                block, blockstart, blockstop, blocklengths, isblockbanded, pointscompatible,
                AbstractProductSpace, MultivariateFun, BivariateSpace,
                @wrapperstructure, @wrapperspaces, @wrapper, @calculus_operator, resizedata!, slnorm, affine_setdiff,
                complexroots
import ApproxFunFourier: PeriodicDomain, PeriodicLine, Circle
import FFTW: plan_r2r!, fftwNumber, REDFT10, REDFT01, REDFT00, RODFT00, R2HC, HC2R,
                r2r!, r2r,  plan_fft, plan_ifft, plan_ifft!, plan_fft!, cFFTWPlan
import Base: values, convert, getindex, setindex!, *, +, -, ==, <, <=, >, |, !, !=, eltype, iterate,
                                >=, /, ^, \, ∪, transpose, size, tail, broadcast, broadcast!, copyto!, copy, to_index, (:),
                                similar, map, vcat, hcat, hvcat, show, summary, stride, sum, cumsum, sign, imag, conj, inv,
                                complex, reverse, exp, sqrt, abs, abs2, sign, issubset, values, in, first, last, rand, intersect, setdiff,
                                isless, union, angle, join, isnan, isapprox, isempty, sort, merge, promote_rule,
                                minimum, maximum, extrema, argmax, argmin, findmax, findmin, isfinite,
                                zeros, zero, one, promote_rule, repeat, length, resize!, isinf,
                                getproperty, findfirst, unsafe_getindex, fld, cld, div, real, imag,
                                @_inline_meta, eachindex, firstindex, lastindex, keys, isreal, OneTo,
                                Array, Vector, Matrix, view, ones, @propagate_inbounds, print_array,
                                split

import ApproxFunOrthogonalPolynomials: Laguerre

import LinearAlgebra: ⋅

export PeriodicLine, chop!, OscLaurent, SumFun, OscConstantSpace, zai, cai, condense, Cauchy, CauchyP, CauchyM, ⋅, fouriertransform, FourierTransform, combine!#, spacescompatible
#include("Domains/Domains.jl")

# Override the evaluation of Piecewise Funs
evaluate(f::Fun{T},x) where {T <: PiecewiseSpace{S}} where S = sum(map(F -> F(x),components(f)))

struct OscLaurent{D<:PeriodicLine,R} <: Space{D,R} # OscLaurent{D<:SPeriodicLine,R}?
    domain::D
    exp::R
    OscLaurent{D,R}(d,ex) where {D,R} = new{D,R}(d,ex)
    OscLaurent{D,R}(d) where {D,R} = new{D,R}(d,0.)
    OscLaurent{D,R}() where {D,R} = new{D,R}(D(),0.,0)
end
OscLaurent(d::PeriodicLine,exp::Float64) = OscLaurent{typeof(d),prectype(d)}(d,exp)
OscLaurent(d::PeriodicLine) = OscLaurent(d,0.)
OscLaurent(α::Float64) = OscLaurent(PeriodicLine{false,Complex{Float64}}(0.,1.),α)
OscLaurent(α::Float64,L::Float64) = OscLaurent(PeriodicLine{false,Complex{Float64}}(0.,L),α)

OscLaurent(d::PeriodicLine,exp::BigFloat) = OscLaurent{typeof(d),prectype(d)}(d,exp)
OscLaurent(α::BigFloat) = OscLaurent(PeriodicLine{false,BigFloat}(0.,1.),α)
OscLaurent(α::BigFloat,L::BigFloat) = OscLaurent(PeriodicLine{false,BigFloat}(0.,L),α)

OscLaurent() = OscLaurent(PeriodicLine())

include("SumFun.jl")

struct OscConstantSpace{D<:PeriodicLine,R} <: Space{D,R} # OscLaurent{D<:SPeriodicLine,R}?
    domain::D
    exp::Float64
    OscConstantSpace{D,R}(d,ex) where {D,R} = new{D,R}(d,ex)
    OscConstantSpace{D,R}(d) where {D,R} = new{D,R}(d,0.)
    OscConstantSpace{D,R}() where {D,R} = new{D,R}(D(),0.,0)
end
OscConstantSpace(d::PeriodicLine,exp::Float64) = OscConstantSpace{typeof(d),complex(prectype(d))}(d,exp)
OscConstantSpace(d::PeriodicLine) = OscConstantSpace(d,0.)
OscConstantSpace(α::Float64) = OscConstantSpace(PeriodicLine{false,Complex{Float64}}(0.,1.),α)
OscConstantSpace(α::Float64,L::Float64) = OscConstantSpace(PeriodicLine{false,Complex{Float64}}(0.,L),α)
OscConstantSpace() = OscConstantSpace(PeriodicLine())

maxspace(a::LaguerreWeight,b::LaguerreWeight) =  spacescompatible(a,b) ? a : PiecewiseSpace(a,b)
maxspace(a::OscLaurent,b::OscConstantSpace) = a
maxspace(b::OscConstantSpace,a::OscLaurent) = a

function getindex(C::ConcreteConversion{A,B,T},k::Integer,j::Integer) where {A<:OscConstantSpace,B<:OscLaurent,T}
    if j == k == 1
        return one(T)
    else
        return zero(T)
    end
end
Base.size(C::ConcreteConversion{A,B,T}) where {A<:OscConstantSpace,B<:OscLaurent,T} = (∞,1)
+(f::Fun{S},g::Fun{T}) where {S<:OscConstantSpace,T<:OscLaurent} = g+f
-(f::Fun{S},g::Fun{T}) where {S<:OscConstantSpace,T<:OscLaurent} = g-f
bandwidths(C::ConcreteConversion{A,B,T}) where {A<:OscConstantSpace,B<:OscLaurent,T} = (0,0)
Conversion(A::OscConstantSpace,B::OscLaurent) = ConcreteConversion(A::OscConstantSpace,B::OscLaurent)

function *(f::Fun{OscLaurent{D,R}},g::Fun{OscLaurent{D,R}})  where {D,R}
    sp = OscLaurent(f.space.domain,f.space.exp+g.space.exp)
    #m = maximum([length(f.coefficients),length(g.coefficients)])
    m = length(f.coefficients) + length(g.coefficients)
    Fun(sp,ApproxFun.transform(sp,values(pad(f,m)).*values(pad(g,m))))
end

function *(f::Fun{OscConstantSpace{D,R}},g::Fun{OscLaurent{D,R}})  where {D,R}
    sp = OscLaurent(f.space.domain,f.space.exp+g.space.exp)
    #m = maximum([length(f.coefficients),length(g.coefficients)])
    #m = length(f.coefficients) + length(g.coefficients)
    Fun(sp,f.coefficients[1]*g.coefficients)
end

function *(f::Fun{OscConstantSpace{D,R}},g::Fun{OscConstantSpace{DD,R}})  where {D,DD,R}
    sp = OscConstantSpace(f.space.domain,f.space.exp+g.space.exp)
    #m = maximum([length(f.coefficients),length(g.coefficients)])
    #m = length(f.coefficients) + length(g.coefficients)
    Fun(sp,f.coefficients[1]*g.coefficients[1])
end


evaluate(f::Fun{OscConstantSpace{D,R}},x::Float64) where {D,R} = f.coefficients[1]*exp(1im*f.space.exp*x)
function *(A::Array{T,2},b::Array{S,1})  where {T<:Fun,S<:Fun}
    c = [];
    for i = 1:size(A)[1]
        sum = A[i,1]*b[1]
        for j = 2:size(A)[2]
            sum = sum + A[i,j]*b[j]
        end
        append!(c,[sum])
    end
    c
end

function *(A::Array{T,2},b::Array{S,1})  where {T<:SumFun,S<:SumFun}
    c = [];
    for i = 1:size(A)[1]
        sum = A[i,1]*b[1]
        for j = 2:size(A)[2]
            sum = sum + A[i,j]*b[j]
        end
        append!(c,[sum])
    end
    c
end

spacescompatible(A::OscConstantSpace{D,R},B::OscLaurent{D,R})  where {D,R} = false
spacescompatible(B::OscLaurent{D,R},A::OscConstantSpace{D,R}) where {D,R} = A.exp ≈ B.exp
spacescompatible(B::OscConstantSpace{D,R},A::OscConstantSpace{D,R}) where {D,R} = A.exp ≈ B.exp
#hasconversion(A::OscConstantSpace,B::OscLaurent) = A.exp == B.exp
#hasconversion(B::OscLaurent,A::OscConstantSpace) = false

## A hack, definitely
## There is an issue because if F = Fun(f,OscLaurent(α)) then F != f, typically
## If checkpoints could depend on space and domain, this could be fixed easily
function Fun(f::Function,sp::OscLaurent{D,R}) where {D,R}
    g = Fun(f,Laurent(sp.domain))
    n = length(g.coefficients)
    n = round(Int64,n/length(g(points(g)[1])))+1
    return Fun(f,sp,n)
end

## Inner product for vectors.  Could extend to matrices using trace.
function ⋅(f::Fun{T},g::Fun{S}) where {S<:ApproxFun.ArraySpace{Q,1},T<:ApproxFun.ArraySpace{J,1}} where {Q,J}
    sum(map(⋅,f,g))
end

function ⋅(f::Fun{T},g::Fun{S},s::Function) where {S<:ApproxFun.ArraySpace{Q,1},T<:ApproxFun.ArraySpace{J,1}} where {Q,J}
    sum(s(transpose(conj(f))*g))[1]
end

function ⋅(f::Fun{T},g::Fun{S},s) where {S,T}
    sum(s(conj(f)*g))
end

function ⋅(f::Fun{T},g::Fun{S}) where {S,T}
    sum(conj(f)*g)
end

function ⋅(f::SumFun ,g::SumFun)
    sum(conj(f)*g)
end

#function ⋅(f::SumFun ,g::SumFun)
#    sum(conj(f)*g)
#end

function ⋅(f::Array{T,1},g::Array{S,1}) where {S<:SumFun,T<:SumFun}
    sum(sum(map(*,conj(f),g)))
end

function ⋅(f::Array{T,1},g::Array{S,1}) where {S<:Fun,T<:Fun}
    sum(map(⋅,f,g))
end

function condense(f::Fun)
    f
end

function condense(f::Array{T,1}) where T<:Fun
    map(condense,f)
end

function condense(f::Fun{T}) where {T <: SumSpace}
    sum(components(f))
end

function condense(f::Fun{T}) where {T <: ArraySpace{S}} where {S<:SumSpace}
    Fun(map( x -> sum(components(x)), Array(f)))
end

function condense(f::Fun{T}) where {T <: ArraySpace{S}} where {S<:Space}
    f
end

function chop(f::Array{T,1}) where T<:Fun
    map(chop,f)
end

function chop!(f::Array{T,1}) where T<:Fun
    map(chop!,f)
end

spacescompatible(a::OscLaurent{D,R},b::OscLaurent{D,R}) where {D,R} = a.exp ≈ b.exp
fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints(::Type{T},n::Integer) where {T<:Number} = convert(T,π)*collect(0:2:2n-2)/n

##FFT That interlaces coefficients
zero_nan(x) = isnan(x) ? zero(typeof(x)) : x

function fsum(x::AbstractVector{T}) where T
    x[1] + sum(x[2:2:end].*((-1)^j for j in 1:length(x[2:2:end]))) + sum(x[3:2:end].*((-1)^j for j in 1:length(x[3:2:end])))
end


struct plan_mfft!
    P
end
plan_mfft!(x::AbstractVector{T}) where T = plan_mfft!(plan_fft!(x))

function *(tr::plan_mfft!,x::AbstractVector{T}) where T
    tr.P*x
    return x
end

struct plan_mifft!
    P
end
plan_mifft!(x::AbstractVector{T}) where T = plan_mifft!(plan_ifft!(x))

function *(tr::plan_mifft!,x::AbstractVector{T}) where T
    return tr.P*x
end

plan_transform!(sp::OscLaurent{D,R},x::AbstractVector{T}) where {D,R,T<:Complex} =
    TransformPlan(sp,plan_fft!(x),Val{true})
plan_itransform!(sp::OscLaurent{D,R},x::AbstractVector{T}) where {D,R,T<:Complex} =
    ITransformPlan(sp,plan_ifft!(x),Val{true})

plan_transform!(sp::OscLaurent{D,R},x::AbstractVector{T}) where {D,R,T} =
    error("In place variants not possible with real data.")
plan_itransform!(sp::OscLaurent{D,R},x::AbstractVector{T}) where {D,R,T} =
    error("In place variants not possible with real data.")

plan_transform(sp::OscLaurent{D,R},x::AbstractVector{T}) where {T<:Complex,D,R} =
    TransformPlan(sp,plan_transform!(sp,x),Val{false})
plan_itransform(sp::OscLaurent{D,R},x::AbstractVector{T}) where {T<:Complex,D,R} =
    ITransformPlan(sp,plan_itransform!(sp,x),Val{false})

function *(P::TransformPlan{T,OscLaurent{D,R},true},vals::AbstractVector{T}) where {T,D,R}
    n = length(vals)
    #vals = [zero_nan(j) for j in vals]
    vals = lmul!(inv(convert(T,n)),P.plan*vals)
    reverseeven!(interlace!(vals,1))
    vals
end

function *(P::ITransformPlan{T,OscLaurent{D,R},true},cfs::AbstractVector{T}) where {T,D,R}
    n = length(cfs)
    reverseeven!(cfs)
    cfs[:]=[cfs[1:2:end];cfs[2:2:end]]  # TODO: deinterlace!
    lmul!(n,cfs)
    P.plan*cfs
end

*(P::TransformPlan{T,OscLaurent{D,R},false},vals::AbstractVector{T}) where {T<:Complex,D,R} =
    P.plan*copy(vals)
*(P::TransformPlan{T,OscLaurent{D,R},false},vals::AbstractVector{T}) where {T,D,R} =
    P.plan*Vector{Complex{T}}(vals)
*(P::ITransformPlan{T,OscLaurent{D,R},false},vals::AbstractVector{T}) where {T<:Complex,D,R} =
    P.plan*copy(vals)
*(P::ITransformPlan{T,OscLaurent{D,R},false},vals::AbstractVector{T}) where {T,D,R} =
    P.plan*Vector{Complex{T}}(vals)

transform(::OscLaurent{D,R},vals,plan) where {D,R} = plan*vals
itransform(::OscLaurent{D,R},cfs,plan) where {D,R} = plan*cfs

transform(sp::OscLaurent{D,R},vals::AbstractVector) where {D,R} = plan_transform(sp,vals)*vals
itransform(sp::OscLaurent{D,R},cfs::AbstractVector) where {D,R} = plan_itransform(sp,cfs)*cfs



include("utils.jl")
include("calculus.jl")
include("operators.jl")
include("cauchy.jl")

function Base.conj(sp::OscLaurent{DD,RR}) where {DD,RR}
    OscLaurent(domain(sp),-sp.exp)
end


function Base.conj(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}  ## Almost the same as Laurent
    ncoefficients(f) == 0 && return f

    cfs = Array{cfstype(f)}(undef, iseven(ncoefficients(f)) ? ncoefficients(f)+1 : ncoefficients(f))
    cfs[1] = conj(f.coefficients[1])
    cfs[ncoefficients(f)] = 0
    for k=2:2:ncoefficients(f)-1
        @inbounds cfs[k] = conj(f.coefficients[k+1])
    end
    for k=3:2:ncoefficients(f)+1
        @inbounds cfs[k] = conj(f.coefficients[k-1])
    end
    Fun(conj(space(f)),cfs)
end

end #module
