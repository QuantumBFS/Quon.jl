using Printf

abstract type AbstractQuonParam end

include("quon_const.jl")

mutable struct QuonParam{T <: Number} <: AbstractQuonParam
    param::T
    isparallel::Bool
end

function Base.show(io::IO, p::QuonParam)
    print(io, p.isparallel ? "∥ " : "⊥ ")
    re = real(p.param)
    ig = imag(p.param)
    if !isapprox(re, 0; atol = quon_atol)
        @printf(io, "%.3f + ", re)
    end
    @printf(io, "%.3f im", ig)
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

function change_direction(p::QuonParam)
    new_param = change_direction(p.param)
    new_isparallel = !p.isparallel
    return QuonParam{typeof(new_param)}(new_param, new_isparallel)
end

function is_singular_change_direction(p::QuonParam)
    (real(p.param) ≈ 0 && rem(imag(p.param), pi, RoundNearest) ≈ 0) && return true
    return false
end

function is_singular_yang_baxter(p1::QuonParam, p2::QuonParam, p3::QuonParam)
    any(is_singular_change_direction(p) for p in (p1, p2, p3)) && return true
    try
        yang_baxter_param(p1, p2, p3)
    catch e
        return true
    end
    return false
end

# function update_yang_baxter_triangle(p1, p2, p3)
#     q1, q2, q3 = copy(p1), copy(p2), copy(p3)
#     if is_param_pi(p1) && !p1.isparallel
#         if is_param_pi(p2) && !p2.isparallel
#             q1.param = zero(p1.param)
#             q2.param = zero(p2.param)
#             if is_param_pi(p3) && !p3.isparallel
#                 q3.param = zero(q3.param)
#             else
#                 q3 = q3 + QuonParam(pi*im, false)
#             end
#         elseif is_param_pi(p3) && !p3.isparallel
#             q1.param = zero(p1.param)
#             q3.param = zero(p3.param)
#             q2 = q2 + QuonParam(pi*im, false)
#         end
#     elseif is_param_pi(p2) && !p2.isparallel && is_param_pi(p3) && !p3.isparallel
#         q2.param = zero(p1.param)
#         q3.param = zero(p3.param)
#         q1 = q1 + QuonParam(pi*im, false)
#     end
#     return q1, q2, q3
# end
# function update_yang_baxter_star(p1, p2, p3)
#     q1, q2, q3 = p1, p2, p3
#     if is_param_pi(p1) && p1.isparallel
#         if is_param_pi(p2) && p2.isparallel
#             q1.param = zero(p1.param)
#             q2.param = zero(p2.param)
#             if is_param_pi(p3) && p3.isparallel
#                 q3.param = zero(q3.param)
#             else
#                 q3 = q3 + QuonParam(pi*im, true)
#             end
#         elseif is_param_pi(p3) && p3.isparallel
#             q1.param = zero(p1.param)
#             q3.param = zero(p3.param)
#             q2 = q2 + QuonParam(pi*im, true)
#         end
#     elseif is_param_pi(p2) && p2.isparallel && is_param_pi(p3) && p3.isparallel
#         q2.param = zero(p1.param)
#         q3.param = zero(p3.param)
#         q1 = q1 + QuonParam(pi*im, true)
#     end
#     return q1, q2, q3
# end

is_parallel(p::QuonParam) = p.isparallel

is_param_pi(p::QuonParam; atol = quon_atol) = (isapprox(real(p.param), 0, atol = atol) && 
    isapprox(rem(imag(p.param)+pi, 2pi, RoundNearest), 0, atol = atol))
is_param_zero(p::QuonParam; atol = quon_atol) = (isapprox(real(p.param), 0, atol = atol) && 
    isapprox(rem(imag(p.param), 2pi, RoundNearest), 0, atol = atol))
is_param_half_pi(p::QuonParam; atol = quon_atol) = (isapprox(real(p.param), 0, atol = atol) && 
    isapprox(rem(imag(p.param), 2pi, RoundNearest), π/2, atol = atol))
is_param_minus_half_pi(p::QuonParam; atol = quon_atol) = (isapprox(real(p.param), 0, atol = atol) && 
    isapprox(rem(imag(p.param), 2pi, RoundNearest), -π/2, atol = atol))

"""
    is_param_para_half_pi(p; atol = atol)

Returns `true` if `p` ≈ `π/2 (∥)` or `-π/2 (⊥)`.
"""
is_param_para_half_pi(p::QuonParam; atol = quon_atol) = (p.isparallel && is_param_half_pi(p; atol = atol)) || 
    (!p.isparallel && is_param_minus_half_pi(p; atol = atol))

"""
    is_param_perp_half_pi(p; atol = atol)

Returns `true` if `p` ≈ `π/2 (⊥)` or `-π/2 (∥)`.
"""
is_param_perp_half_pi(p::QuonParam; atol = quon_atol) = (!p.isparallel && is_param_half_pi(p; atol = atol)) || 
    (p.isparallel && is_param_minus_half_pi(p; atol = atol))

function yang_baxter_param(p1::QuonParam, p2::QuonParam, p3::QuonParam)
    # triangle => star
    p1.isparallel && (p1 = change_direction(p1))
    !p2.isparallel && (p2 = change_direction(p2))
    p3.isparallel && (p3 = change_direction(p3))
    
    q1_param, q2_param, q3_param = yang_baxter_param(p1.param, p2.param, p3.param)
    q1, q2, q3 = QuonParam(q1_param, true), QuonParam(q2_param, false), QuonParam(q3_param, true)
    return q1, q2, q3
end
function yang_baxter_param_inv(p1::QuonParam, p2::QuonParam, p3::QuonParam)
    # star => triangle
    !p1.isparallel && (p1 = change_direction(p1))
    p2.isparallel && (p2 = change_direction(p2))
    !p3.isparallel && (p3 = change_direction(p3))

    q1_param, q2_param, q3_param = yang_baxter_param_inv(p1.param, p2.param, p3.param)
    return QuonParam(q1_param, false), QuonParam(q2_param, true), QuonParam(q3_param, false)
end