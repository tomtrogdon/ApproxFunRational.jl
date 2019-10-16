
## Differentiation
# Only k = 1 for now...
Derivative(S::OscLaurent{DD,RR},k::Integer) where {DD,RR} = ConcreteDerivative(S,k)
bandwidths(D::ConcreteDerivative{OscLaurent{DD,RR}})  where {DD,RR} = (3,3)
rangespace(D::ConcreteDerivative{S}) where {S<:OscLaurent}=D.space
function mob_derivative_getindex(d::PeriodicLine,m,j::Integer,k::Integer)
    # j is the row index, i.e. (j,k) not (k,j) as other ApproxFun implementations
    if j == 1
        0.0im
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

## Multiplication, not same as Laurent
Multiplication(f::Fun{OscLaurent{DD,RR}},sp::OscLaurent{DD,RR}) where {DD,RR} = ConcreteMultiplication(cfstype(f),f,sp)

function osc_laurent_getindex(v::AbstractVector{T},k::Integer,j::Integer) where T
    # switch to double-infinite indices
    negative = v[3:2:end];
    nonnegative = v[[1;2:2:end]]
    l=k;p=j;
    k=iseven(k) ? -k÷2 : (k-1)÷2
    j=iseven(j) ? -j÷2 : (j-1)÷2
    c = sum(negative.*((-1)^j for j in 1:length(negative))) + sum(nonnegative.*((-1)^(j-1) for j in 1:length(nonnegative)))
    b = (l <= length(v)) ? (-1)^j*v[l] : zero(T)
    if j-k == 0
        nonnegative[1] - c - b
    elseif l == 1 && p > 1
        zero(T) - b 
    elseif 0<k-j≤length(negative)
        negative[k-j] - b
    elseif 0<j-k≤length(nonnegative)-1
        nonnegative[j-k+1] - b
    else
        zero(T) - b
    end
end

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
    #laurent_getindex(T.f.coefficients[3:2:end],T.f.coefficients[[1;2:2:end]],k,j)
    osc_laurent_getindex(T.f.coefficients,k,j)
end

function bandwidths(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}}) where {DD,RR}
    bbi = blockbandwidths(T)
    (2bbi[1]+1,2bbi[2]+1)
end

function blockbandwidths(T::ConcreteMultiplication{OscLaurent{DD,RR},OscLaurent{DD,RR}}) where {DD,RR}
    m = ncoefficients(T.f)÷2
    (m,m)
end
