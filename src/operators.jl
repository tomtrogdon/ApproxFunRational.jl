
## Differentiation
# Only k = 1 for now...
Derivative(S::OscLaurent{DD,RR},k::Integer) where {DD,RR} = ConcreteDerivative(S,k)
bandwidths(D::ConcreteDerivative{OscLaurent{DD,RR}})  where {DD,RR} = (3,3)
rangespace(D::ConcreteDerivative{S}) where {S<:OscLaurent}=D.space
function mob_derivative_getindex(d::PeriodicLine,m,j::Integer,k::Integer)
    # j is the row index, i.e. (j,k) not (k,j) as other ApproxFun implementations
    if j == 1 && 2 <= k <= 3
        -(-1)^k*0.5im/d.L
    elseif iseven(j)
        l = j÷2
        if k == j
            -1.0im*l/d.L
        elseif k == j+2
            -1.0im*(l+1)/2*1/d.L
        elseif k == j-2 && j > 2
            -1.0im*(l-1)/2*1/d.L
        else
            0.0im
        end
    else
        l = j÷2
        if k == j
            1.0im*l/d.L
        elseif k == j+2
            1.0im*(l+1)/2*1/d.L
        elseif k == j-2 && j > 3
            1.0im*(l-1)/2*1/d.L
        else
            0.0im
        end
    end
end

getindex(D::ConcreteDerivative{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} =
convert(T,(k == j ? domainspace(D).exp*1im : 0.) + mob_derivative_getindex(domain(D),D.order,k,j))

## Multiplication, same as Laurent
Multiplication(f::Fun{OscLaurent{DD,RR}},sp::OscLaurent{DD,RR}) where {DD,RR} = ConcreteMultiplication(cfstype(f),f,sp)

function laurent_getindex(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer) where T
    # switch to double-infinite indices
    k=iseven(k) ? -k÷2 : (k-1)÷2
    j=iseven(j) ? -j÷2 : (j-1)÷2

    if 0<k-j≤length(negative)
        negative[k-j]
    elseif 0≤j-k≤length(nonnegative)-1
        nonnegative[j-k+1]
    else
        zero(T)
    end
end

rangespace(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}}) where {DD,RR} = OscLaurent(domainspace(T).domain,T.space.exp + space(T.f).exp)
function getindex(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}},k::Integer,j::Integer) where {DD,RR}
    isempty(T.f.coefficients) && return zero(eltype(T))
    laurent_getindex(T.f.coefficients[3:2:end],T.f.coefficients[[1;2:2:end]],k,j)

end

function bandwidths(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}}) where {DD,RR}
    bbi = blockbandwidths(T)
    (2bbi[1]+1,2bbi[2]+1)
end

function blockbandwidths(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}}) where {DD,RR}
    m = ncoefficients(T.f)÷2
    (m,m)
end

function fouriertransform(f::Fun{T}) where T <: OscLaurent
    α = f.space.exp
    L = f.space.domain.L
    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
    c1 = f.coefficients[3:2:end]
    c2 = f.coefficients[2:2:end]
    c1 = c1.*[4*π*L*(-1)^(i+1) for i=1:length(c1)]
    c2 = c2.*[4*π*L*(-1)^(i+1) for i=1:length(c2)]
    Fun(LaguerreWeight(0.0,0.5,sp1),c1) + Fun(LaguerreWeight(0.0,0.5,sp2),c2)
end

#Op class = FourierOperator
#Op = FourierTransform
#ConcOp = ConcreteFourierTransform
#WrappOp = FourierTransformWrapper

abstract type FourierOperator{S,OT,T} <: Operator{T} end

abstract type FourierTransform{S,OT,T} <: FourierOperator{S,OT,T} end

struct ConcreteFourierTransform{S<:Space,OT,T} <: FourierTransform{S,OT,T}
    space::S
    sign::OT
end
struct FourierTransformWrapper{BT<:Operator,S<:Space,OT,T} <: FourierTransform{S,OT,T}
    op::BT
    sign::OT
end

@wrapper FourierTransformWrapper

ConcreteFourierTransform(sp::Space,k) = ConcreteFourierTransform{typeof(sp),typeof(k),prectype(sp)}(sp,k)

FourierTransform(sp::UnsetSpace,k) = ConcreteFourierTransform(sp,k)
ConcreteFourierTransform(sp::Space) = ConcreteFourierTransform(sp,1.0)

FourierTransform(sp::OscLaurent,k) = ConcreteFourierTransform(sp,k)
# not needed yet
function Base.convert(::Type{Operator{T}},D::ConcreteFourierTransform) where T
    if T==eltype(D)
        D
    else
        ConcreteFourierTransform{typeof(D.space),T}(D.space)
    end
end

function Base.convert(::Type{Operator{T}},D::FourierTransformWrapper) where T
    if T==eltype(D)
        D
    else
        # work around typeinfernece bug
        op=convert(Operator{T},D.op)
        FourierTransformWrapper{typeof(op),typeof(domainspace(op)),T}(op)
    end
end

ApproxFunBase.domain(D::ConcreteFourierTransform) = domain(D.space)
ApproxFunBase.domainspace(D::ConcreteFourierTransform) = D.space
Base.getindex(::ConcreteFourierTransform{UnsetSpace,T},k::Integer,j::Integer) where {OT,T} =
    error("Spaces cannot be inferred for operator")
ApproxFunBase.rangespace(D::ConcreteFourierTransform{UnsetSpace,OT,T}) where {OT,T} = UnsetSpace()
ApproxFunBase.promotedomainspace(D::FourierTransform,sp::UnsetSpace) = D
function ApproxFunBase.promotedomainspace(D::FourierTransform,sp::Space)
    if isambiguous(domain(sp))
        Fouriertransform(typeof(sp)(domain(D)))
    else
        FourierTransform(sp)
    end
end

choosedomainspace(M::FourierOperator{UnsetSpace},sp::Space) =
    iswrapper(M) ? choosedomainspace(M.op,sp) : sp

function rangespace(D::ConcreteFourierTransform{S,OT,T}) where {S<:OscLaurent,OT,T}
    α = D.space.exp
    L = D.space.domain.L
    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
    SumSpace(LaguerreWeight(0.0,0.5,sp1),LaguerreWeight(0.0,0.5,sp2))
end

function osclaurent_ft_getindex(k::Integer,j::Integer,L::Float64,T::Type)
    if j != k - 1
        return zero(T)
    else
        # switch to double-infinite indices
        # k=iseven(k) ? -k÷2 : (k-1)÷2
        j = iseven(j) ? -j÷2 : (j-1)÷2
        return 4*π*L*(-1)^(j+1)
    end
end

getindex(D::ConcreteFourierTransform{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} = osclaurent_ft_getindex(k,j,D.space.domain.L,T)
bandwidths(D::ConcreteFourierTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = (0,1)
