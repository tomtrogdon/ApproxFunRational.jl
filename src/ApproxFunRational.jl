module ApproxFunRational
using Base, LinearAlgebra, Reexport, AbstractFFTs, FFTW, InfiniteArrays, FillArrays, FastTransforms, IntervalSets,
            DomainSets, SpecialFunctions

@reexport using ApproxFunBase

import ApproxFunFourier: PeriodicDomain, PeriodicLine, Circle

export OscLaurent, ic_s!, c_s!, mobiusdiff, zai, cai#, spacescompatible
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

spacescompatible(a::OscLaurent{D,R},b::OscLaurent{D,R}) where {D,R} = a.exp == b.exp

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints(::Type{T},n::Integer) where {T<:Number} = convert(T,π)*collect(0:2:2n-2)/n

##FFT That interlaces coefficients
zero_nan(x) = isnan(x) ? zero(typeof(x)) : x

function fsum(x::AbstractVector{T}) where T
    x[1] + sum(x[2:2:end].*((-1)^j for j in 1:length(x[2:2:end]))) + sum(x[3:2:end].*((-1)^j for j in 1:length(x[3:2:end])))
end

c_s!(x) = x[1] = fsum(x) # normalize
ic_s!(x) = x[1] = x[1] - (fsum(x)-x[1])

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
    ic_s(x)
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
include("Operators.jl")




end #module
