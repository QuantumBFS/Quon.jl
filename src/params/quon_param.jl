using Printf

abstract type AbstractQuonParam end

mutable struct QuonParam{T} <: AbstractQuonParam
    param::T
    isparallel::Bool
end

function Base.show(io::IO, p::QuonParam{T}) where T
    print(io, p.isparallel ? "∥ " : "⊥ ")
    if p isa Complex
        re = real(p.param)
        ig = imag(p.param)
        if !isapprox(re, 0; atol = quon_atol)
            @printf(io, "%.3f + ", re)
        end
        @printf(io, "%.3f⋅im", ig)
    elseif p isa QuonConst
        p === Zero && print(io, "0")
        p === Pi && print(io, "π⋅im")
        p === HalfPi && print(io, "π/2⋅im")
        p === NegHalfPi && print(io, "-π/2⋅im")
        p === InfZero && print(io, "CD(0)")
        p === InfPi && print(io, "CD(π⋅im)")
    else
        print(io, p)
    end
end

Base.copy(p::QuonParam{T}) where T = QuonParam{T}(p.param, p.isparallel)

function Base.:(+)(p1::QuonParam, p2::QuonParam)
    if p1.isparallel == p2.isparallel
        new_param = p1.param + p2.param
        new_ispara = p1.isparallel
    elseif !isinf(change_direction(p1.param))
        new_param = change_direction(p1.param) + p2.param
        new_ispara = p2.isparallel
    elseif !isinf(change_direction(p2.param))
        new_param = p1.param + change_direction(p2.param)
        new_ispara = p1.isparallel
    else
        return !p1.isparallel ? p1 : p2
    end
    return QuonParam{typeof(new_param)}(new_param, new_ispara)
end

function change_direction(p::QuonParam{T}) where {T <: QuonComplex}
    new_param = change_direction(p.param)
    new_isparallel = !p.isparallel
    return QuonParam{T}(new_param, new_isparallel)
end

function is_singular_yang_baxter(p1::QuonParam{T}, p2::QuonParam{T}, p3::QuonParam{T}) where {T <:QuonComplex}
    try
        yang_baxter_param(p1, p2, p3)
    catch e
        return true
    end
    return false
end

is_parallel(p::QuonParam) = p.isparallel

is_pi(p::Complex; atol = quon_atol) = (isapprox(real(p), 0, atol = atol) && 
    isapprox(rem(imag(p)+pi, 2pi, RoundNearest), 0, atol = atol))
is_zero(p::Complex; atol = quon_atol) = (isapprox(real(p), 0, atol = atol) && 
    isapprox(rem(imag(p), 2pi, RoundNearest), 0, atol = atol))
is_half_pi(p::Complex; atol = quon_atol) = (isapprox(real(p), 0, atol = atol) && 
    isapprox(rem(imag(p), 2pi, RoundNearest), π/2, atol = atol))
is_neg_half_pi(p::Complex; atol = quon_atol) = (isapprox(real(p), 0, atol = atol) && 
    isapprox(rem(imag(p), 2pi, RoundNearest), -π/2, atol = atol))
is_inf_zero(p::Complex) = (real(p) == -Inf)
is_inf_pi(p::Complex) = (real(p) == Inf)

is_para_pi(p::QuonParam) = (is_parallel(p) && is_pi(p)) || (!is_parallel(p) && is_inf_pi(p))
is_perp_pi(p::QuonParam) = (!is_parallel(p) && is_pi(p)) || (is_parallel(p) && is_inf_pi(p))
is_para_zero(p::QuonParam) = is_parallel(p) && is_zero(p)
is_perp_zero(p::QuonParam) = !is_parallel(p) && is_zero(p)


"""
    is_para_half_pi(p)

Returns `true` if `p` ≈ `π/2 (∥)` or `-π/2 (⊥)`.
"""
is_para_half_pi(p::QuonParam) = (is_parallel(p) && is_half_pi(p)) || 
    (!is_parallel(p) && is_neg_half_pi(p))

"""
    is_perp_half_pi(p)

Returns `true` if `p` ≈ `π/2 (⊥)` or `-π/2 (∥)`.
"""
is_perp_half_pi(p::QuonParam) = (!is_parallel(p) && is_half_pi(p)) || 
    (is_parallel(p) && is_neg_half_pi(p))

function yang_baxter_param(p1::QuonParam{T}, p2::QuonParam{T}, p3::QuonParam{T}) where {T <: QuonComplex}
    # triangle => star
    is_parallel(p1) && (p1 = change_direction(p1))
    !is_parallel(p2) && (p2 = change_direction(p2))
    is_parallel(p3) && (p3 = change_direction(p3))
    
    q1_param, q2_param, q3_param = yang_baxter_param(p1.param, p2.param, p3.param)
    q1_param = to_quon_const(q1_param)
    q2_param = to_quon_const(q2_param)
    q3_param = to_quon_const(q3_param)
    q1, q2, q3 = QuonParam{T}(q1_param, true), QuonParam{T}(q2_param, false), QuonParam{T}(q3_param, true)
    return q1, q2, q3
end
function yang_baxter_param_inv(p1::QuonParam{T}, p2::QuonParam{T}, p3::QuonParam{T}) where {T <: QuonComplex}
    # star => triangle
    !is_parallel(p1) && (p1 = change_direction(p1))
    is_parallel(p2) && (p2 = change_direction(p2))
    !is_parallel(p3) && (p3 = change_direction(p3))

    q1_param, q2_param, q3_param = yang_baxter_param_inv(p1.param, p2.param, p3.param)
    q1_param = to_quon_const(q1_param)
    q2_param = to_quon_const(q2_param)
    q3_param = to_quon_const(q3_param)
    return QuonParam{T}(q1_param, false), QuonParam{T}(q2_param, true), QuonParam{T}(q3_param, false)
end