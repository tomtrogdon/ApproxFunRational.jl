
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

# old code for testing
#function fouriertransform(f::Fun{T}) where T <: OscLaurent
#    α = f.space.exp
#    L = f.space.domain.L
#    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
#    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
#    c1 = f.coefficients[3:2:end]
#    c2 = f.coefficients[2:2:end]
#    c1 = c1.*[4*π*L*(-1)^(i+1) for i=1:length(c1)]
#    c2 = c2.*[4*π*L*(-1)^(i+1) for i=1:length(c2)]
#    Fun(LaguerreWeight(0.0,0.5,sp1),c1) + Fun(LaguerreWeight(0.0,0.5,sp2),c2)
#end

#Op class = FourierOperator
#Op = FourierTransform
#ConcOp = ConcreteFourierTransform
#WrappOp = FourierTransformWrapper


### Fourier transform

abstract type FourierOperator{S,OT,T} <: Operator{T} end

abstract type FourierTransform{S,OT,T} <: FourierOperator{S,OT,T} end

for Op in (:FourierTransform,:δFourierTransform)
    ConcOp = Meta.parse("Concrete"*string(Op))
    OpWrap = Meta.parse(string(OP)*"Wrapper")

    struct $ConcOp{S<:Space,OT,T} <: FourierTransform{S,OT,T}
        space::S
        sign::OT
    end

    struct $OpWrap{BT<:Operator,S<:Space,OT,T} <: FourierTransform{S,OT,T}
        op::BT
        sign::OT
    end

    @wrapper $OpWrap

    $ConcOp(sp::Space,k) = $ConcOp{typeof(sp),typeof(k),prectype(sp)}(sp,k)
    $ConcOp(sp::Space) = $ConcOp(sp,1.0)


    # not needed yet
    function Base.convert(::Type{Operator{T}},D::$ConcOp) where T
        if T==eltype(D)
            D
        else
            $ConcOp{typeof(D.space),T}(D.space)
        end
    end

    function Base.convert(::Type{Operator{T}},D::$WrapOp) where T
        if T==eltype(D)
            D
        else
            # work around typeinfernece bug
            op=convert(Operator{T},D.op)
            $WrapOp{typeof(op),typeof(domainspace(op)),T}(op)
        end
    end

    ApproxFunBase.domain(D::$ConcOp) = domain(D.space)
    ApproxFunBase.domainspace(D::$ConcOp) = D.space
    Base.getindex(::$ConcOp{UnsetSpace,T},k::Integer,j::Integer) where {OT,T} =
        error("Spaces cannot be inferred for operator")
        ApproxFunBase.rangespace(D::$ConcOp{UnsetSpace,OT,T}) where {OT,T} = UnsetSpace()
ApproxFunBase.promotedomainspace(D::FourierTransform,sp::UnsetSpace) = D
function ApproxFunBase.promotedomainspace(D::FourierTransform,sp::Space)
    if isambiguous(domain(sp))
        FourierTransform(typeof(sp)(domain(D)))
    else
        FourierTransform(sp)
    end
end

choosedomainspace(M::FourierOperator{UnsetSpace},sp::Space) =
    iswrapper(M) ? choosedomainspace(M.op,sp) : sp

function rangespace(D::ConcreteFourierTransform{S,OT,T}) where {S<:OscLaurent,OT,T}
    α = D.space.exp/abs(D.sign)
    L = D.space.domain.L*abs(D.sign)
    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
    if D.sign > 0.
        SumSpace(LaguerreWeight(0.0,0.5,sp1),LaguerreWeight(0.0,0.5,sp2))
    else
        SumSpace(LaguerreWeight(0.0,0.5,sp2),LaguerreWeight(0.0,0.5,sp1))
    end
end

FourierTransform(sp::UnsetSpace,k) = ConcreteFourierTransform(sp,k)
FourierTransform(sp::OscLaurent,k) = ConcreteFourierTransform(sp,k)
FourierTransform(sp::OscLaurent) = ConcreteFourierTransform(sp,1.0)
FourierTransform() = ConcreteFourierTransform(UnsetSpace(),1.0)

function osclaurent_ft_getindex(k::Integer,j::Integer,L::Float64,T::Type)
    if j != k + 1
        return zero(T)
    else
        # switch to double-infinite indices
        # k=iseven(k) ? -k÷2 : (k-1)÷2
        j = iseven(j) ? -j÷2 : (j-1)÷2
        return 4*π*L*(-1)^(j+1)
    end
end

getindex(D::ConcreteFourierTransform{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} = osclaurent_ft_getindex(k,j,D.space.domain.L,T)
bandwidths(D::ConcreteFourierTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = 0,1
#blockbandwidths(D::ConcreteFourierTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = ((0,∞),(0,1))
### Cauchy transform

abstract type SingularIntegralOperator{S,OT,T} <: Operator{T} end

abstract type CauchyTransform{S,OT,T} <: SingularIntegralOperator{S,OT,T} end

struct ConcreteFourierTransform{S<:Space,T} <: FourierTransform{S,OT,T}
    space::S
end
struct FourierTransformWrapper{BT<:Operator,S<:Space,T} <: FourierTransform{S,OT,T}
    op::BT
end

@wrapper CauchyTransformWrapper

ConcreteCauchyTransform(sp::Space) = ConcreteCauchyTransform{typeof(sp),prectype(sp)}(sp)

CauchyTransform(sp::UnsetSpace) = ConcreteCauchyTransform(sp)
ConcreteCauchyTransform(sp::Space) = ConcreteFourierTransform(sp)

CauchyTransform(sp::OscLaurent) = ConcreteCauchyTransform(sp)
CauchyTransform() = ConcreteCauchyTransform(UnsetSpace())
# not needed yet
function Base.convert(::Type{Operator{T}},D::ConcreteFourierTransform) where T
    if T==eltype(D)
        D
    else
        ConcreteFourierTransform{typeof(D.space),T}(D.space)
    end
end

function Base.convert(::Type{Operator{T}},D::CauchyTransformWrapper) where T
    if T==eltype(D)
        D
    else
        # work around typeinfernece bug
        op=convert(Operator{T},D.op)
        CauchyTransformWrapper{typeof(op),typeof(domainspace(op)),T}(op)
    end
end

ApproxFunBase.domain(D::ConcreteCauchyTransform) = domain(D.space)
ApproxFunBase.domainspace(D::ConcreteCauchyTransform) = D.space
Base.getindex(::ConcreteCauchyTransform{UnsetSpace,T},k::Integer,j::Integer) where {OT,T} =
    error("Spaces cannot be inferred for operator")
ApproxFunBase.rangespace(D::ConcreteCauchyTransform{UnsetSpace,T}) where {OT,T} = UnsetSpace()
ApproxFunBase.promotedomainspace(D::CauchyTransform,sp::UnsetSpace) = D
function ApproxFunBase.promotedomainspace(D::CauchyTransform,sp::Space)
    if isambiguous(domain(sp))
        CauchyTransform(typeof(sp)(domain(D)))
    else
        CauchyTransform(sp)
    end
end

choosedomainspace(M::CauchyOperator{UnsetSpace},sp::Space) =
    iswrapper(M) ? choosedomainspace(M.op,sp) : sp

function rangespace(D::ConcreteCauchyTransform{S,T}) where {S<:OscLaurent,T}
    if D.space.exp != 0.0
        dom = D.space
        SumSpace(D.space,OscLaurent(dom,0.0))
    else
        D.space
    end
end

function osclaurent_cauchy_getindex(k::Integer,j::Integer,L::Float64,T::Type)
    if j != k + 1
        return zero(T)
    else
        # switch to double-infinite indices
        # k=iseven(k) ? -k÷2 : (k-1)÷2
        j = iseven(j) ? -j÷2 : (j-1)÷2
        return 4*π*L*(-1)^(j+1)
    end
end

getindex(D::ConcreteCauchyTransform{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} = osclaurent_ft_getindex(k,j,D.space.domain.L,T)
bandwidths(D::ConcreteCauchyTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = (1,1)
