using Printf

abstract type AbstractQuonParam end

mutable struct QuonParam{T} <: AbstractQuonParam
    param::T
    isparallel::Bool
end

function Base.show(io::IO, p::QuonParam{T}) where T
    print(io, p.isparallel ? "∥ " : "⊥ ")
    if p.param isa Complex
        re = real(p.param)
        ig = imag(p.param)
        if !isapprox(re, 0; atol = quon_atol)
            @printf(io, "%.3f + ", re)
        end
        @printf(io, "%.3f⋅im", ig)
    elseif p.param isa QuonConst
        p.param === Zero && print(io, "0")
        p.param === Pi && print(io, "π⋅im")
        p.param === HalfPi && print(io, "π/2⋅im")
        p.param === NegHalfPi && print(io, "-π/2⋅im")
        p.param === InfZero && print(io, "CD(0)")
        p.param === InfPi && print(io, "CD(π⋅im)")
    else
        print(io, p.param)
    end
end

Base.copy(p::QuonParam{T}) where T = QuonParam{T}(p.param, p.isparallel)

function add_with_dir(p1::QuonParam{T}, p2::QuonParam{T}, is_para::Bool) where T
    is_parallel(p1) == is_para || (p1 = change_direction(p1))
    is_parallel(p2) == is_para || (p2 = change_direction(p2))
    return QuonParam{T}(p1.param + p2.param, is_para)
end

add_parallel(p1::QuonParam{T}, p2::QuonParam{T}) where T = add_with_dir(p1, p2, true)
add_orthorgonal(p1::QuonParam{T}, p2::QuonParam{T}) where T = add_with_dir(p1, p2, false)

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

is_para_pi(p::QuonParam) = (is_parallel(p) && is_pi(p.param)) || 
    (!is_parallel(p) && is_inf_pi(p.param))
is_orth_pi(p::QuonParam) = (!is_parallel(p) && is_pi(p.param)) || 
    (is_parallel(p) && is_inf_pi(p.param))
is_para_zero(p::QuonParam) = (is_parallel(p) && is_zero(p.param)) || 
    (!is_parallel(p) && is_inf_zero(p.param))
is_orth_zero(p::QuonParam) = (!is_parallel(p) && is_zero(p.param)) || 
    (is_parallel(p) && is_inf_zero(p.param))


"""
    is_para_half_pi(p)

Returns `true` if `p` ≈ `π/2 (∥)` or `-π/2 (⊥)`.
"""
is_para_half_pi(p::QuonParam) = (is_parallel(p) && is_half_pi(p.param)) || 
    (!is_parallel(p) && is_neg_half_pi(p.param))

"""
    is_orth_half_pi(p)

Returns `true` if `p` ≈ `π/2 (⊥)` or `-π/2 (∥)`.
"""
is_orth_half_pi(p::QuonParam) = (!is_parallel(p) && is_half_pi(p.param)) || 
    (is_parallel(p) && is_neg_half_pi(p.param))

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