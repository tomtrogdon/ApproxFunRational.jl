module ApproxFunRational
using Base, ApproxFun, ApproxFunBase, ApproxFunFourier, Reexport, FFTW, LinearAlgebra#, Reexport, AbstractFFTs, FFTW, InfiniteArrays, FillArrays, FastTransforms, IntervalSets,
            #DomainSets, SpecialFunctions

@reexport using ApproxFunBase

#import ApproxFunBase: Fun, AbstractTransformPlan, TransformPlan, ITransformPlan, ConcreteDerivative, ConcreteMultiplication, prectype, spacescompatible

import ApproxFunBase: normalize!, flipsign, FiniteRange, Fun, MatrixFun, UnsetSpace, VFun, RowVector,
                UnivariateSpace, AmbiguousSpace, SumSpace, SubSpace, WeightSpace, NoSpace, Space,
                HeavisideSpace, PointSpace,
                IntervalOrSegment, RaggedMatrix, AlmostBandedMatrix,
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
                Integral,
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



export PeriodicLine, OscLaurent, zai, cai, CauchyPNO, CauchyP, CauchyM#, spacescompatible
#include("Domains/Domains.jl")

struct OscLaurent{D<:PeriodicLine,R} <: Space{D,R} # OscLaurent{D<:SPeriodicLine,R}?
    domain::D
    exp::Float64
    OscLaurent{D,R}(d,ex) where {D,R} = new{D,R}(d,ex)
    OscLaurent{D,R}(d) where {D,R} = new{D,R}(d,0.)
    OscLaurent{D,R}() where {D,R} = new{D,R}(D(),0.,0)
end
OscLaurent(d::PeriodicLine,exp::Float64) = OscLaurent{typeof(d),complex(prectype(d))}(d,exp)
OscLaurent(d::PeriodicLine) = OscLaurent(d,0.)
OscLaurent() = OscLaurent(PeriodicLine())

# need to overload for adaptivity
function Fun(f::Function,sp::OscLaurent{D,R}) where {D,R}
    g = Fun(f,Laurent(sp.domain))
    return Fun(sp,g.coefficients)
end

spacescompatible(a::OscLaurent{D,R},b::OscLaurent{D,R}) where {D,R} = a.exp == b.exp

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints(::Type{T},n::Integer) where {T<:Number} = convert(T,π)*collect(0:2:2n-2)/n

##FFT That interlaces coefficients
zero_nan(x) = isnan(x) ? zero(typeof(x)) : x

function fsum(x::AbstractVector{T}) where T
    x[1] + sum(x[2:2:end].*((-1)^j for j in 1:length(x[2:2:end]))) + sum(x[3:2:end].*((-1)^j for j in 1:length(x[3:2:end])))
end

c_s!(x) = 1# do nothing %x[1] = fsum(x) # normalize
ic_s!(x) = 1# do nothing x[1] = x[1] - (fsum(x)-x[1])

struct plan_mfft!
    P
end
plan_mfft!(x::AbstractVector{T}) where T = plan_mfft!(plan_fft!(x))

function *(tr::plan_mfft!,x::AbstractVector{T}) where T
    tr.P*x
    c_s!(x)
    return x
end

struct plan_mifft!
    P
end
plan_mifft!(x::AbstractVector{T}) where T = plan_mifft!(plan_ifft!(x))

function *(tr::plan_mifft!,x::AbstractVector{T}) where T
    ic_s!(x)
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
    c_s!(vals)
    vals
end

function *(P::ITransformPlan{T,OscLaurent{D,R},true},cfs::AbstractVector{T}) where {T,D,R}
    n = length(cfs)
    ic_s!(cfs)
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
        cfs[k] = conj(f.coefficients[k+1])
    end
    for k=3:2:ncoefficients(f)+1
        cfs[k] = conj(f.coefficients[k-1])
    end
    Fun(conj(space(f)),cfs)
end

end #module
