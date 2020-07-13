function cai(f,c)
    return x -> isnan(f(x)) ? c : complex(f(x))
end

function zai(f)
    return x -> isnan(f(x)) ? 0.0im : complex(f(x))
end

function horner(c::AbstractVector,kr::AbstractRange{Int},x)
    T = promote_type(eltype(c),eltype(x))
    if isempty(c)
        return zero(x)
    end
    ret = zero(T)
    @inbounds for k in reverse(kr)
        ret = muladd(x,ret,c[k])
    end
    ret
end

function horner(c::AbstractVector,kr::AbstractRange{Int},x::AbstractVector)
    n,T = length(x),promote_type(eltype(c),eltype(x))
    if isempty(c)
        return zero(x)
    end

    ret = zeros(T,n)
    @inbounds for k in reverse(kr)
        ck = c[k]
        @simd for i = 1:n
            ret[i] = muladd(x[i],ret[i],ck)
        end
    end
    ret
end

horner(c::AbstractVector,x) = horner(c,1:length(c),x)
horner(c::AbstractVector,x::AbstractArray) = horner(c,1:length(c),x)
horner(c::AbstractVector,kr::AbstractRange{Int},x::AbstractArray) = reshape(horner(c,kr,vec(x)),size(x))

function evaluate(f::AbstractVector,S::OscLaurent{D,R},z) where {D,R}
    #z = tocanonical(domain(S),z)
    zz = mappoint(domain(S),Circle(),z)
    invz = 1/zz
    (horner(f,1:2:length(f),zz) + horner(f,2:2:length(f),invz).*invz)*exp.(1im*z*S.exp)
end

function mobiusdiff(v::AbstractVector{T}) where T<:Number #need to fix
    n = length(v)
    if n == 1
        w = zeros(T,1)
    else
        w = Array{T}(undef, n+2)
        w[1] = zero(T)
        n = length(v)
        k = 1 #just for clarity
        @inbounds w[3] = 1.0im*k*v[3] + 0.5im*(k+1)*v[5]
        @inbounds w[2] = -1.0im*k*v[2] - 0.5im*(k+1)*v[4]
        l = (isodd(n) ? n÷2 : n÷2-1)
        for k=2:l-1
            @inbounds w[2k] = -1.0im*k*v[2k] - 0.5im*(k+1)*v[2k+2] - 0.5im*(k-1)*v[2k-2]
            @inbounds w[2k+1] = 1.0im*k*v[2k+1] + 0.5im*(k+1)*v[2k+3] + 0.5im*(k-1)*v[2k-1]
        end
        #### Deal with large n case
        if iseven(n) # if n is even, then l = n/2 -1, n = 2l + 2
            # positive modes
            @inbounds w[2l] = -1.0im*l*v[2l] - 0.5im*(l+1)*v[2l+2] - 0.5im*(l-1)*v[2l-2]
            @inbounds w[2l+2] = -1.0im*(l+1)*v[2l+2] - 0.5im*(l)*v[2l]
            @inbounds w[2l+4] = -0.5im*(l+1)*v[2l+2] # index -(n/2+2)
            # negative modes
            @inbounds w[2l+1] = 1.0im*l*v[2l+1] + 0.5im*(l-1)*v[2l-1]
            @inbounds w[2l+3] = 0.5im*(l)*v[2l+1] # index n/2
        else # if n is odd, then things are symmetric
            @inbounds w[2l+1] = 1.0im*l*v[2l+1] + 0.5im*(l-1)*v[2l-1]
            @inbounds w[2l+3] = 0.5im*(l)*v[2l+1]
            @inbounds w[2l] = -1.0im*l*v[2l] - 0.5im*(l-1)*v[2l-2]
            @inbounds w[2l+2] = -0.5im*(l)*v[2l]
        end
        w[1] = -fsum(w) # needed because fsum(w) = 0 is needed to enforce decay at infinity
        w
    end
end

function mob_conj(x::AbstractVector{T}) where T <: Number
    y = x
    if isodd(length(x))
        y[2:2:end] = conj(x[3:2:end])
        y[3:2:end] = conj(x[2:2:end])
    else  # This is not the exact conjugate!
        y[2:2:end-1] = conj(x[3:2:end-1])
        y[3:2:end-1] = conj(x[2:2:end-1])
    end
end
