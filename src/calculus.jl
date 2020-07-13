differentiate(f::Fun{OscLaurent{DD,RR}}) where {DD,RR} = (f.space.exp ≈ 0.0) ?  Fun(f.space,1/(domain(f.space).L)*mobiusdiff(f.coefficients)) : Fun(f.space,vcat(1.0im*f.space.exp*f.coefficients,[0.0,0.0]) + 1/(domain(f.space).L)*mobiusdiff(f.coefficients))


function Lag(n::Int64,x::Float64) # evaluate Laguerre polynomials at x, α = 1
  c = zeros(Float64,n+2)
  c[2] = exp.(-x/2)
  k = 0
  for i = 3:n+2 #k = 1 is the zeroth order, gives (k-1)th order
    c[i] = (2*k+2-x)*c[i-1]/(k+1) - c[i-2]
    k = k + 1
  end
  return c
end

Base.sum(f::Fun{T}) where T <: ArraySpace = map(sum,f)


# Old version, still here for testing
function sum2(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}
    if space(f).exp ≈ 0.
        data = f.coefficients[2:2:end]
        data = data.*(j*(-1)^j for j in 1:length(data))
        s = sum(data)
        data = f.coefficients[3:2:end]
        data = data.*(j*(-1)^j for j in 1:length(data))
        return -2*pi*(s + sum(data))*domain(f).L
    else
        m = convert(Int64,floor(length(coefficients(f))/2))
        ex = space(f).exp
        β = domain(f).L
        z = 2*β*abs(ex)
        c = Lag(m-1,z)[2:end]
        #cc = v.coefs
        #if abs(v.exp) == 0
        #  c = c.*(N(2*m)[2:length(c)+1])
        #end # if v.exp > 0, only negative terms contribute
        if ex > 0.
          data = f.coefficients[2:2:end]
        elseif ex < 0.
            data = f.coefficients[3:2:end]
        end
        data = data.*((-1)^j for j in 1:length(data))
        return -4*pi*β*sum(data.*c[1:length(data)])
    end
end

function Base.sum(f::Fun{OscLaurent{DD,RR}}) where {DD,RR <: BigFloat}
    if space(f).exp ≈ 0
        s = BigFloat("0")*1im
        o = 1
        j = 0
        for i = 2:2:length(f.coefficients)
            o *= -1
            j += 1
            @inbounds s += f.coefficients[i]*j*o
        end
        o = 1
        j = 0
        for i = 3:2:length(f.coefficients)
            o *= -1
            j += 1
            @inbounds s += f.coefficients[i]*j*o
        end
        return -2*big(pi)*s*BigFloat(domain(f).L)
    else
        #m = convert(Int64,floor(length(coefficients(f))/2))
        ex = space(f).exp
        β = BigFloat(domain(f).L)
        z = 2*β*abs(ex)
        c0 = 0
        c1 = exp(-z/2)
        #c = Lag(m-1,z)[2:end]
        #cc = v.coefs
        #if abs(v.exp) == 0
        #  c = c.*(N(2*m)[2:length(c)+1])
        #end # if v.exp > 0, only negative terms contribute
        if real(ex) > 0
          #data = f.coefficients[2:2:end]
          j = 2
        else
            #data = f.coefficients[3:2:end]
            j = 3
        end
        s = BigFloat("0")*1im
        o = 1
        k = 0
        for i = j:2:length(f.coefficients)
            o *= -1
            c2 = (2*k +2 -z)*c1/(k+1) - c0
            @inbounds s += o*c1*f.coefficients[i]
            c0 = c1
            c1 = c2
            k += 1
        end

        #data = data.*((-1)^j for j in 1:length(data))
        return -4*big(pi)*β*s
    end
end

function Base.sum(f::Fun{OscLaurent{DD,RR}}) where {DD,RR}
    if space(f).exp ≈ 0.
        s = 0.0im
        o = 1
        j = 0
        for i = 2:2:length(f.coefficients)
            o *= -1
            j += 1
            @inbounds s += f.coefficients[i]*j*o
        end
        o = 1
        j = 0
        for i = 3:2:length(f.coefficients)
            o *= -1
            j += 1
            @inbounds s += f.coefficients[i]*j*o
        end
        return -2*pi*s*domain(f).L
    else
        #m = convert(Int64,floor(length(coefficients(f))/2))
        ex = space(f).exp
        β = domain(f).L
        z = 2*β*abs(ex)
        c0 = 0.
        c1 = exp(-z/2)
        #c = Lag(m-1,z)[2:end]
        #cc = v.coefs
        #if abs(v.exp) == 0
        #  c = c.*(N(2*m)[2:length(c)+1])
        #end # if v.exp > 0, only negative terms contribute
        if real(ex) > 0.
          #data = f.coefficients[2:2:end]
          j = 2
        else
            #data = f.coefficients[3:2:end]
            j = 3
        end
        s = 0.0im
        o = 1
        k = 0
        for i = j:2:length(f.coefficients)
            o *= -1
            c2 = (2*k +2 -z)*c1/(k+1) - c0
            @inbounds s += o*c1*f.coefficients[i]
            c0 = c1
            c1 = c2
            k += 1
        end

        #data = data.*((-1)^j for j in 1:length(data))
        return -4*pi*β*s
    end
end
