mutable struct SumFun <: Function# OscLaurent{D<:SPeriodicLine,R}?
    funs::Vector{Fun}
end

evaluate(f::SumFun,x) = sum(map(y -> y(x),f.funs))

SumFun(f::Fun{T}) where T = SumFun([f])
SumFun(f::Fun{T}) where T<: SumSpace = SumFun(components(f))
SumFun(f::Fun{T}) where T<: PiecewiseSpace = SumFun(components(f))

conj(s::SumFun) = SumFun(map(conj,s.funs))
conj(s::Array{T}) where T<:SumFun = map(conj,s)

#sum(s::SumFun) = sum(map(sum,s.funs))
sum(s::SumFun) = s.funs == [] ? 0.0 : sum(map(sum,s.funs))
sum(s::Array{T}) where T<:SumFun = map(sum,s)

function combine!(s::SumFun)
    s.funs = combine!(s.funs)
    s
end
combine!(s::Vector{SumFun}) = map(combine!,s)
combine!(s::Array{SumFun}) = map(combine!,s)
components(f::Fun) = f

function combine!(f::Vector{T}) where T <: Fun
    f = vcat(map(components,f))
    j = 0
    @inbounds while j < length(f)
        j += 1
        if (length(f[j].coefficients) <= 4 && sum(abs.(f[j].coefficients)) < 1e-15 && length(f) > 2) || maximum(abs.(f[j].coefficients)) < 1e-15
            deleteat!(f,j)
            j -= 1
        end
    end
    i = 0
    @inbounds while i < length(f)
        i += 1
        j = i
        while j < length(f)
            j += 1
            if spacescompatible(f[i],f[j])
                f[i] += f[j]
                deleteat!(f,j)
                j -= 1
            end
        end
    end
    f
end

+(f::SumFun,g::SumFun) = SumFun(vcat(g.funs,f.funs))
-(f::SumFun,g::SumFun) = SumFun(vcat(-g.funs,f.funs))

+(f::SumFun,g::Fun) = f + SumFun(g)
+(g::Fun,f::SumFun) = f + SumFun(g)

-(f::SumFun,g::Fun) = f - SumFun(g)
-(g::Fun,f::SumFun) = -f + SumFun(g)

-(f::SumFun) = SumFun(map(y -> -y,f.funs))

*(f::Fun,g::SumFun) = SumFun(map(y -> f*y,g.funs))
*(f::SumFun,g::Fun) = SumFun(map(y -> g*y,f.funs))
*(f::SumFun,g::Number) = SumFun(map(y -> g*y,f.funs))
*(f::Number,g::SumFun) = SumFun(map(y -> f*y,g.funs))
*(f::SumFun,g::Array) = map(x -> f*x,g)

*(f::SumFun,g::SumFun) = combine!(+(map(y -> y*g,f.funs)...))

*(c::Float64,g::SumFun) = SumFun(map(y -> c*y,g.funs))
*(c::Complex{Float64},g::SumFun) = SumFun(map(y -> c*y,g.funs))

(f::SumFun)(x) = evaluate(f,x)
+(f::SumFun) = f

(f::Array)(x) = map(y -> evaluate(y,x),f)
chop!(f::SumFun) = SumFun(map(chop!,f.funs))
chop!(f::Array) = map(chop!,f)

chop!(f::SumFun,tol::Float64) = SumFun(map(x -> chop!(x,tol),f.funs))
chop!(f::Array,tol::Float64) = map(x -> chop!(x,tol),f)
copy(f::SumFun) = SumFun(copy(f.funs))
