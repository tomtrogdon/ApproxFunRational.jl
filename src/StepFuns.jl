struct OscLaurentStep{D<:PeriodicLine,R} <: Space{D,R} # OscLaurent{D<:SPeriodicLine,R}?
    domain::D
    exp::Float64
    OscLaurentStep{D,R}(d,ex) where {D,R} = new{D,R}(d,ex)
    OscLaurentStep{D,R}(d) where {D,R} = new{D,R}(d,0.)
    OscLaurentStep{D,R}() where {D,R} = new{D,R}(D(),0.,0)
end
OscLaurentStep(d::PeriodicLine,exp::Float64) = OscLaurentStep{typeof(d),complex(prectype(d))}(d,exp)
OscLaurentStep(d::PeriodicLine) = OscLaurentStep(d,0.)
OscLaurentStep(α::Float64) = OscLaurentStep(PeriodicLine{false,Float64}(0.,1.),α)
OscLaurentStep(α::Float64,L::Float64) = OscLaurentStep(PeriodicLine{false,Float64}(0.,L),α)
OscLaurentStep() = OscLaurentStep(PeriodicLine())

blocklengths(C::OscLaurentStep) = [1]
block(C::OscLaurentStep,k) = Block(1)
