function LagSeries(j::Int64,z::Float64,x::AbstractVector{T}) where T # need to investigate stability
  c = Lag(j,z) # c[2] gives L_0^(1), c[i] gives L_{i-2}^(1)
  v0 = x
  out = zeros(Complex{Float64},length(x),j)
  out[1:end,1] = v0*c[2]
  for i = 2:j
    out[1:end,i] = (1.0 .+ x).*out[1:end,i-1] .+ x*(c[i+1]-c[i])
  end

  return out
end

# old version, still here for testing
function LagSeries_old(z::Float64,x::AbstractVector{T},cfs::AbstractVector{S}) where {T,S} # need to investigate stability
  j = length(cfs)
  c = Lag(j,z) # c[2] gives L_0^(1), c[i] gives L_{i-2}^(1)
  v0 = x
  #out = zeros(Complex{Float64},length(x))
  out = v0*c[2]
  pm = -1
  sum = (pm*cfs[1])*out
  for i = 2:j
    pm = pm*(-1)
    @inbounds out = (1.0 .+ x).*out .+ x*(c[i+1]-c[i])
    @inbounds sum = sum + (pm*cfs[i])*out
  end
  return sum
end

function LagSeries(z::Float64,x::AbstractVector{T},cfs::AbstractVector{S}) where {T,S} # need to investigate stability
  j = length(cfs)
  #c = Lag(j,z) # c[2] gives L_0^(1), c[i] gives L_{i-2}^(1)
  #v0 = x
  #out = zeros(Complex{Float64},length(x))
  c0 = 0.
  c1 = exp(-z/2) # c[2]
  out = x*c1
  k = 0
  pm = -1
  sum = (pm*cfs[1])*out
  for i = 2:j
    pm = pm*(-1)
    c2 = (2*k+2-z)*c1/(k+1) - c0 # c2 = c[i + 1]
    #out = (1.0 .+ x).*out .+ x*(c2-c1)
    out .*= (1.0 .+ x)
    out .+= x*(c2-c1)
    #@inbounds sum = sum + (pm*cfs[i])*out
    @inbounds axpy!(pm*cfs[i],out,sum)
    k += 1
    c0 = c1
    c1 = c2
  end
  return sum
end

function Res(j::Int64,α::Float64,β::Float64,z::AbstractVector{T}) where T<:Number
  x = -2im*sign(j)*β./(z.+1im*sign(j)*β)
  y = -2*sign(j)*α*β
  return -LagSeries(abs(j),y,x)
end

function Res(j::Integer,α::Float64,β::Float64,z::AbstractVector{T},cfs::AbstractVector{S}) where {T<:Number,S<:Number}
  x = -2im*sign(j)*β./(z.+1im*sign(j)*β)
  y = -2*sign(j)*α*β
  return -LagSeries(y,x,cfs)
end

function CauchyPNO(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}# returns the non-oscillatory portion, assumes decay
    α = space(f).exp
    sp = OscLaurent(domain(space(f)),0.)
    if α == 0.
      b = copy(f.coefficients)
      b[[1;2:2:end]] = zeros(typeof(f.coefficients[1]),length(b[[1;2:2:end]])) # zero positive
      b[1] = -fsum(b)
      return Fun(sp,b)
    elseif α < 0.
      #m = length(f.coefficients[3:2:end])
      #c = Res(m,α,domain(space(f)).L,points(f))
      return Fun(sp,ApproxFun.transform(sp,-Res(1,α,domain(space(f)).L,points(f),f.coefficients[3:2:end])))
      #return Fun(sp,ApproxFun.transform(sp,-c*([(-1)^i for i = 1:m].*f.coefficients[3:2:end])))
    else
      #m = length(f.coefficients[2:2:end])
      #c = Res(-m,α,domain(space(f)).L,points(f))
      return Fun(sp,ApproxFun.transform(sp,Res(-1,α,domain(space(f)).L,points(f),f.coefficients[2:2:end])))
      #return Fun(sp,ApproxFun.transform(sp,c*([(-1)^i for i = 1:m].*f.coefficients[2:2:end])))
    end
end

function CauchyP(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}
  α = space(f).exp
  if α <= 0.
    return CauchyPNO(f)
  else
    return CauchyPNO(f) + f
  end
end

function CauchyM(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}
  α = space(f).exp
  if α >= 0.
    return CauchyMNO(f)
  else
    return CauchyMNO(f) - f
  end
end


function CauchyMNO(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}# returns the non-oscillatory portion
    α = space(f).exp
    sp = OscLaurent(domain(space(f)),0.)
    if α == 0.
      b = -copy(f.coefficients)
      b[[1;3:2:end]] = zeros(typeof(f.coefficients[1]),length(b[[1;3:2:end]])) # zero positive
      b[1] = -fsum(b)
      return Fun(sp,b)
    elseif α < 0.
      #m = length(f.coefficients[3:2:end])
      #c = Res(m,α,domain(space(f)).L,points(f))
      #return Fun(sp,ApproxFun.transform(sp,-c*([(-1)^i for i = 1:m].*f.coefficients[3:2:end])))
      return Fun(sp,ApproxFun.transform(sp,-Res(1,α,domain(space(f)).L,points(f),f.coefficients[3:2:end])))
    else
      #m = length(f.coefficients[2:2:end])
      #c = Res(-m,α,domain(space(f)).L,points(f))
      #return Fun(sp,ApproxFun.transform(sp,c*([(-1)^i for i = 1:m].*f.coefficients[2:2:end])))
      return Fun(sp,ApproxFun.transform(sp,Res(-1,α,domain(space(f)).L,points(f),f.coefficients[2:2:end])))
    end
end

for cauchy in (:CauchyP,:CauchyM)
  @eval begin
    function $cauchy(f::Fun{T}) where {T<:SumSpace}
      sum(map(x -> $cauchy(x),components(f)))
    end

    function $cauchy(f::Fun{T}) where {T<:PiecewiseSpace}
      sum(map(x -> $cauchy(x),components(f)))
    end

    function $cauchy(f::Fun{T}) where {T<:ApproxFun.ArraySpace{J,1}} where {J}
      Fun(map($cauchy,Array(f)))
    end
  end
end

abstract type SIOperator{S,OT,T} <: Operator{T} end

macro SI_operator(Op)
    ConcOp = Meta.parse("Concrete"*string(Op))
    OpWrap = Meta.parse(string(Op)*"Wrapper")
    return esc(quote
        abstract type $Op{SSS,OT,TTT} <: SIOperator{SSS,OT,TTT} end

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
        $ConcOp(sp::Space) = $ConcOp(sp,1)

        $Op(sp::UnsetSpace,k) = $ConcOp(sp,k)
        $Op(sp::UnsetSpace,k::Integer) = $ConcOp(sp,k)

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


#Cauchy(σ,f) = σ == 1 ? CauchyP(f) : CauchyM(f)
@SI_operator(CauchyOperator)

function rangespace(D::ConcreteCauchyOperator{S,OT,T}) where {S<:OscLaurent,OT,T}
    α = D.space.exp
    L = D.space.domain.L
    return SumSpace(OscLaurent(α,L),OscLaurent(0.0,L))
end
bandwidths(D::ConcreteCauchyOperator{OscLaurent{DD,RR},OT,T}) where {DD,RR,OT,T} = ∞,∞
getindex(D::ConcreteCauchyOperator{OscLaurent{DD,RR},OT,T},k::Integer,j::Integer) where {DD,RR,OT,T} = one(T)

##TODO: Finish implementation of CauchyOperator



abstract type AbstractCauchyOperator end

struct Cauchy <: AbstractCauchyOperator
  pm::Integer
end

*(C::AbstractCauchyOperator,F::Fun) = (C.pm == 1) ? CauchyP(F) : CauchyM(F)
