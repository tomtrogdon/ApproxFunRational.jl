
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


### Fourier transform

abstract type FourierOperator{S,OT,T} <: Operator{T} end

macro fourier_operator(Op)
    ConcOp = Meta.parse("Concrete"*string(Op))
    OpWrap = Meta.parse(string(Op)*"Wrapper")
    return esc(quote
        abstract type $Op{SSS,OT,TTT} <: FourierOperator{SSS,OT,TTT} end

        struct $ConcOp{S<:Space,OT,T} <: $Op{S,OT,T}
            space::S
            sign::OT
        end

        struct $OpWrap{BT<:Operator,S<:Space,OT,T} <: $Op{S,OT,T}
            op::BT
            sign::OT
        end

        @wrapper $OpWrap

        $ConcOp(sp::Space,k) = $ConcOp{typeof(sp),typeof(k),prectype(sp)}(sp,k)
        $ConcOp(sp::Space) = $ConcOp(sp,1.0)

        $Op(sp::UnsetSpace,k) = $ConcOp(sp,k)
        $Op(sp::UnsetSpace,k::Real) = $ConcOp(sp,k)

        # not needed yet
        function Base.convert(::Type{Operator{T}},D::$ConcOp) where T
            if T==eltype(D)
            D
            else
                $ConcOp{typeof(D.space),typeof(D.sign),T}(D.space,D.sign)
            end
        end

        function Base.convert(::Type{Operator{T}},D::$OpWrap) where T
            if T==eltype(D)
            D
            else
                # work around typeinfernece bug
                op=convert(Operator{T},D.op)
                $OpWrap{typeof(op),typeof(domainspace(op)),T}(op)
            end
        end

        ApproxFunBase.domain(D::$ConcOp) = domain(D.space)
        ApproxFunBase.domainspace(D::$ConcOp) = D.space
        Base.getindex(::$ConcOp{UnsetSpace,T},k::Integer,j::Integer) where {OT,T} =
            error("Spaces cannot be inferred for operator")
        ApproxFunBase.rangespace(D::$ConcOp{UnsetSpace,OT,T}) where {OT,T} = UnsetSpace()

        ApproxFunBase.promotedomainspace(D::$Op,sp::UnsetSpace) = D
        function ApproxFunBase.promotedomainspace(D::$Op,sp::Space)
            if isambiguous(domain(sp))
                $Op(typeof(sp)(domain(D)),D.sign)
            else
                $Op(sp,D.sign)
            end
        end
    end)
end
hasconversion(a::LaguerreWeight,b::DiracSpace) = false
hasconversion(b::DiracSpace,a::LaguerreWeight) = false
#maxspace(a::SumSpace,b::DiracSpace) = SumSpace(a,b)
#maxspace(b::DiracSpace,a::SumSpace) = SumSpace(b,a)

@fourier_operator(δFourierTransform)
@fourier_operator(SFourierTransform)

function FourierTransform(k::Float64)
    SFourierTransform(k) + δFourierTransform(k)
end

SFourierTransform(sp::OscLaurent,k) = ConcreteSFourierTransform(sp,k)
SFourierTransform(k::Float64) = ConcreteSFourierTransform(UnsetSpace(),k)
SFourierTransform() = ConcreteSFourierTransform(UnsetSpace(),1.0)

δFourierTransform(sp::OscLaurent,k) = ConcreteδFourierTransform(sp,k)
δFourierTransform(k::Float64) = ConcreteδFourierTransform(UnsetSpace(),k)
δFourierTransform() = ConcreteδFourierTransform(UnsetSpace(),1.0)

choosedomainspace(M::FourierOperator{UnsetSpace},sp::Space) =
    iswrapper(M) ? choosedomainspace(M.op,sp) : sp

function rangespace(D::ConcreteSFourierTransform{S,OT,T}) where {S<:OscLaurent,OT,T}
    α = D.space.exp/abs(D.sign)
    L = D.space.domain.L*abs(D.sign)
    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
    sp0 = DiracSpace(D.space.exp/D.sign)
    if D.sign > 0.
        SumSpace(sp0,LaguerreWeight(0.0,0.5,sp2),LaguerreWeight(0.0,0.5,sp1))
    else
        SumSpace(sp0,LaguerreWeight(0.0,0.5,sp1),LaguerreWeight(0.0,0.5,sp2))
    end
end

function rangespace(D::ConcreteδFourierTransform{S,OT,T}) where {S<:OscLaurent,OT,T}
    α = D.space.exp/abs(D.sign)
    L = D.space.domain.L*abs(D.sign)
    sp1 = Laguerre(1.0,Ray(α,0.0,1/(2.0*L),true))
    sp2 = Laguerre(1.0,Ray(α,π,1/(2.0*L),true))
    sp0 = DiracSpace(D.space.exp/D.sign)
    if D.sign > 0.
        SumSpace(sp0,LaguerreWeight(0.0,0.5,sp2),LaguerreWeight(0.0,0.5,sp1))
    else
        SumSpace(sp0,LaguerreWeight(0.0,0.5,sp1),LaguerreWeight(0.0,0.5,sp2))
    end
end

#TODO: Check this

function rangespace(D::ConcreteδFourierTransform{S,OT,T}) where {S<:LaguerreWeight,OT,T}
    @assert S.α == 0.0 && S.L == 0.5
    sp = D.space
    α = sp.space.domain.center*D.sign
    L = 1/(D.sign*2*sp.space.domain.L)
    OscLaurent(α,L)
end

function rangespace(D::ConcreteSFourierTransform{S,OT,T}) where {S<:LaguerreWeight,OT,T}
    @assert S.α == 0.0 && S.L == 0.5
    sp = D.space
    α = sp.space.domain.center*D.sign
    L = 1/(D.sign*2*sp.space.domain.L)
    OscLaurent(α,L)
end

## TODO: Implement OscConstantSpace to deal with the fourier transform of delta functions

#rangespace(D::ConcreteδFourierTransform{S,OT,T}) where {S<:OscLaurent,OT,T} = DiracSpace(D.space.exp/D.sign)
#Base.size(A::ConcreteδFourierTransform{S,OT,T},k::Integer)  where {S<:OscLaurent,OT,T} = k==1 ? 1 : dimension(domainspace(A))
bandwidths(D::ConcreteδFourierTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = 0,∞

function osclaurent_ft_getindex(k::Integer,j::Integer,L::Float64,T::Type)
    if k == 1 || j != k
        return zero(T)
    else
        # switch to double-infinite indices
        # k=iseven(k) ? -k÷2 : (k-1)÷2
        j = iseven(j) ? -j÷2 : (j-1)÷2
        return -4*π*L*(-1)^(j)
    end
end

getindex(D::ConcreteSFourierTransform{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} = osclaurent_ft_getindex(k,j,D.space.domain.L,T)
bandwidths(D::ConcreteSFourierTransform{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = 0,0
maxspace(a::LaguerreWeight,b::LaguerreWeight) =  spacescompatible(a,b) ? a : PiecewiseSpace(a,b)

function getindex(D::ConcreteδFourierTransform{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T}
    if k != 1
        return zero(T)
    else
        j = iseven(j) ? -j÷2 : (j-1)÷2
        return (-1)^j
    end
 end
