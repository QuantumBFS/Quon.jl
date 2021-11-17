@enum QuonConst begin
    Zero
    Pi
    HalfPi
    NegHalfPi
    InfZero
    InfPi
end

const QuonComplex = Union{QuonConst, ComplexF64}

function Base.Complex(p::QuonConst)
    p === Zero && return 0.0im
    p === Pi && return π*im
    p === HalfPi && return π/2*im
    p === NegHalfPi && return -π/2*im
    @warn "Representing infinity greatness with a Complex number will cause numerical unstable"
    p === InfZero && return -Inf + 0.0im
    p === InfPi && return Inf + 0.0im
end

function Base.:(-)(p::QuonConst)
    p === HalfPi && return NegHalfPi
    p === NegHalfPi && return HalfPi
    p === InfZero && return InfPi
    p === InfPi && return InfZero
    return p
end

function Base.:(+)(p1::QuonConst, p2::QuonConst)
    p1 === Zero && return p2
    p2 === Zero && return p1
    if p1 === Pi
        p2 === Pi && return Zero
        p2 === HalfPi && return NegHalfPi
        p2 === NegHalfPi && return HalfPi
    elseif p1 === HalfPi
        p2 === Pi && return NegHalfPi
        p2 === HalfPi && return Pi
        p2 === NegHalfPi && return Zero
    elseif p1 === NegHalfPi
        p2 === Pi && return HalfPi
        p2 === HalfPi && return Zero
        p2 === NegHalfPi && return Pi
    end
    
    (p1 in (InfZero, InfPi)) && return p1
    return p2
end
function Base.:(+)(p1::Complex, p2::QuonConst)
    p = p1
    p2 === Zero && return to_quon_const(p)
    p2 === Pi && return to_quon_const(p + π*im)
    p2 === HalfPi && return to_quon_const(p + π/2*im)
    p2 === NegHalfPi && return to_quon_const(p - π/2*im)
    return p2
end
Base.:(+)(p1::QuonConst, p2::Complex) = p2 + p1

Base.zero(::QuonConst) = Zero
Base.zero(::Type{QuonConst}) = Zero
Base.zero(::Type{QuonComplex}) = Zero
function Base.imag(p::QuonConst)
    p === Pi && return pi
    p === HalfPi && return pi/2
    p === NegHalfPi && return -pi/2
    return 0.0
end
Base.real(p::QuonConst) = real(Complex(p))

function change_direction(p::QuonConst)
    p === Zero && return InfZero
    p === InfZero && return Zero
    p === Pi && return InfPi
    p === InfPi && return Pi
    p === HalfPi && return NegHalfPi
    p === NegHalfPi && return HalfPi
end


is_pi(p::QuonConst) = (p === Pi)
is_zero(p::QuonConst) = (p === Zero)
is_half_pi(p::QuonConst) = (p === HalfPi)
is_neg_half_pi(p::QuonConst) = (p === NegHalfPi)
is_inf_pi(p::QuonConst) = (p === InfPi)
is_inf_zero(p::QuonConst) = (p === InfZero)

to_quon_const(p::QuonConst) = p
function to_quon_const(p::Complex)
    is_pi(p) && return Pi
    is_zero(p) && return Zero
    is_half_pi(p) && return HalfPi
    is_neg_half_pi(p) && return NegHalfPi
    is_inf_pi(p) && return InfPi
    is_inf_zero(p) && return InfZero

    return p
end

# one const
function yang_baxter_param(α::QuonConst, β::QuonComplex, γ::QuonComplex)
    if α === Zero
        return to_quon_const(β), to_quon_const(γ), Zero
    elseif α === Pi
        return to_quon_const(-β), γ+Pi, Zero
    elseif α === InfZero
        return Zero, InfZero, change_direction(change_direction(β) + γ)
    elseif α === InfPi
        return Pi, InfZero, change_direction(-change_direction(β) + γ)
    end
    return yang_baxter_param(Complex(α), β, γ)
end

function yang_baxter_param(α::QuonComplex, β::QuonConst, γ::QuonComplex)
    if β === Zero
        return Zero, α+γ, Zero
    elseif β === Pi
        return Pi, -α+γ, Zero
    elseif β === InfZero
        return change_direction(α), InfZero, change_direction(γ)
    elseif β === InfPi
        return -change_direction(α), InfZero, -change_direction(γ)
    end
    return yang_baxter_param(α, Complex(β), γ)
end

function yang_baxter_param(α::QuonComplex, β::QuonComplex, γ::QuonConst)
    c, b, a = yang_baxter_param(γ, β, α)
    return a, b, c
end

# two consts
function yang_baxter_param(α::QuonConst, β::QuonConst, γ::QuonComplex)
    if α === Zero
        return to_quon_const(β), to_quon_const(γ), Zero
    elseif α === Pi
        return -β, γ+Pi, Zero
    elseif α === InfZero
        return Zero, InfZero, change_direction(change_direction(β) + γ)
    elseif α === InfPi
        return Pi, InfZero, change_direction(-change_direction(β) + γ)
    end
    if β === Zero
        return Zero, α+γ, Zero
    elseif β === Pi
        return Pi, -α+γ, Zero
    elseif β === InfZero
        return change_direction(α), InfZero, change_direction(γ)
    elseif β === InfPi
        return -change_direction(α), InfZero, -change_direction(γ)
    end
    return yang_baxter_param(Complex(α), Complex(β), γ)
end

function yang_baxter_param(α::QuonComplex, β::QuonConst, γ::QuonConst)
    c, b, a = yang_baxter_param(γ, β, α)
    return a, b, c
end

function yang_baxter_param(α::QuonConst, β::QuonComplex, γ::QuonConst)
    if α === Zero
        return to_quon_const(β), γ, Zero
    elseif α === Pi
        return to_quon_const(-β), γ+Pi, Zero
    elseif α === InfZero
        return Zero, InfZero, change_direction(change_direction(β) + γ)
    elseif α === InfPi
        return Pi, InfZero, change_direction(-change_direction(β) + γ)
    end
    if γ === Zero
        return Zero, α, to_quon_const(β)
    elseif γ === Pi
        return Zero, α+Pi, to_quon_const(-β)
    elseif γ === InfZero
        return change_direction(change_direction(β) + α), InfZero, Zero
    elseif γ === InfPi
        return change_direction(-change_direction(β) + α), InfZero, Pi
    end
    return yang_baxter_param(Complex(α), β, Complex(γ))
end

# three consts
function yang_baxter_param(α::QuonConst, β::QuonConst, γ::QuonConst)
    if (α, β, γ) ⊆ (HalfPi, NegHalfPi)
        if α === HalfPi && β === HalfPi && γ === HalfPi
            return HalfPi, HalfPi, HalfPi
        elseif α === HalfPi && β === HalfPi && γ === NegHalfPi
            return NegHalfPi, HalfPi, HalfPi
        elseif α === HalfPi && β === NegHalfPi && γ === HalfPi
            return NegHalfPi, HalfPi, NegHalfPi
        elseif α === HalfPi && β === NegHalfPi && γ === NegHalfPi
            return NegHalfPi, NegHalfPi, HalfPi
        elseif α === NegHalfPi && β === HalfPi && γ === HalfPi
            return HalfPi, HalfPi, NegHalfPi
        elseif α === NegHalfPi && β === HalfPi && γ === NegHalfPi
            return HalfPi, NegHalfPi, HalfPi
        elseif α === NegHalfPi && β === NegHalfPi && γ === HalfPi
            return HalfPi, NegHalfPi, NegHalfPi
        elseif α === NegHalfPi && β === NegHalfPi && γ === NegHalfPi
            return NegHalfPi, NegHalfPi, NegHalfPi
        end
    end
    if α === Zero
        return β, γ, Zero
    elseif α === Pi
        return -β, γ+Pi, Zero
    elseif α === InfZero
        return Zero, InfZero, change_direction(change_direction(β) + γ)
    elseif α === InfPi
        return Pi, InfZero, change_direction(-change_direction(β) + γ)
    end
    if β === Zero
        return Zero, α+γ, Zero
    elseif β === Pi
        return Pi, -α+γ, Zero
    elseif β === InfZero
        return change_direction(α), InfZero, change_direction(γ)
    elseif β === InfPi
        return -change_direction(α), InfZero, -change_direction(γ)
    end
    if γ === Zero
        return Zero, α, β
    elseif γ === Pi
        return Zero, α+Pi, -β
    elseif γ === InfZero
        return change_direction(change_direction(β) + α), InfZero, Zero
    elseif γ === InfPi
        return change_direction(-change_direction(β) + α), InfZero, Pi
    end
end
