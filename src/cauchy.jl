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

function LagSeries(z::Float64,x::AbstractVector{T},cfs::AbstractVector{S}) where {T,S} # need to investigate stability
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

function Res(j::Int64,α::Float64,β::Float64,z::AbstractVector{T}) where T<:Number
  x = -2im*sign(j)*β./(z.+1im*sign(j)*β)
  y = -2*sign(j)*α*β
  return -LagSeries(abs(j),y,x)
end

function Res(α::Float64,β::Float64,z::AbstractVector{T},cfs::AbstractVector{S}) where {T<:Number,S<:Number}
  j = length(cfs)
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
      return Fun(sp,ApproxFun.transform(sp,-Res(α,domain(space(f)).L,points(f),f.coefficients[3:2:end])))
      #return Fun(sp,ApproxFun.transform(sp,-c*([(-1)^i for i = 1:m].*f.coefficients[3:2:end])))
    else
      #m = length(f.coefficients[2:2:end])
      #c = Res(-m,α,domain(space(f)).L,points(f))
      return Fun(sp,ApproxFun.transform(sp,Res(α,domain(space(f)).L,points(f),f.coefficients[2:2:end])))
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
      b[[1;3:2:end]] = zeros(typeof(f.coefficients[1]),length(b[[1;2:2:end]])) # zero positive
      b[1] = -fsum(b)
      return Fun(sp,b)
    elseif α < 0.
      #m = length(f.coefficients[3:2:end])
      #c = Res(m,α,domain(space(f)).L,points(f))
      #return Fun(sp,ApproxFun.transform(sp,-c*([(-1)^i for i = 1:m].*f.coefficients[3:2:end])))
      return Fun(sp,ApproxFun.transform(sp,-Res(α,domain(space(f)).L,points(f),f.coefficients[3:2:end])))
    else
      #m = length(f.coefficients[2:2:end])
      #c = Res(-m,α,domain(space(f)).L,points(f))
      #return Fun(sp,ApproxFun.transform(sp,c*([(-1)^i for i = 1:m].*f.coefficients[2:2:end])))
      return Fun(sp,ApproxFun.transform(sp,Res(α,domain(space(f)).L,points(f),f.coefficients[2:2:end])))
    end
end
