mutable struct Phase{T <: Number}
    param::T
    isparallel::Bool
end

Base.copy(p::Phase{T}) where T = Phase{T}(p.param, p.isparallel)
function Base.:(+)(p1::Phase, p2::Phase)
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
    return Phase{typeof(new_param)}(new_param, new_ispara)
end

function change_direction(p::Phase)
    new_param = change_direction(p.param)
    new_isparallel = !p.isparallel
    return Phase{typeof(new_param)}(new_param, new_isparallel)
end

function is_singular_change_direction(p::Phase)
    (real(p.param) ≈ 0 && rem(imag(p.param), pi, RoundNearest) ≈ 0) && return true
    return false
end

function is_singular_yang_baxter(p1::Phase, p2::Phase, p3::Phase)
    any(is_singular_change_direction(p) for p in (p1, p2, p3)) && return true
    try
        yang_baxter_param(p1, p2, p3)
    catch e
        return true
    end
    return false
end

function update_yang_baxter_triangle(p1, p2, p3)
    q1, q2, q3 = copy(p1), copy(p2), copy(p3)
    if is_phase_pi(p1) && !p1.isparallel
        if is_phase_pi(p2) && !p2.isparallel
            q1.param = zero(p1.param)
            q2.param = zero(p2.param)
            if is_phase_pi(p3) && !p3.isparallel
                q3.param = zero(q3.param)
            else
                q3 = q3 + Phase(pi*im, false)
            end
        elseif is_phase_pi(p3) && !p3.isparallel
            q1.param = zero(p1.param)
            q3.param = zero(p3.param)
            q2 = q2 + Phase(pi*im, false)
        end
    elseif is_phase_pi(p2) && !p2.isparallel && is_phase_pi(p3) && !p3.isparallel
        q2.param = zero(p1.param)
        q3.param = zero(p3.param)
        q1 = q1 + Phase(pi*im, false)
    end
    return q1, q2, q3
end
function update_yang_baxter_star(p1, p2, p3)
    q1, q2, q3 = p1, p2, p3
    if is_phase_pi(p1) && p1.isparallel
        if is_phase_pi(p2) && p2.isparallel
            q1.param = zero(p1.param)
            q2.param = zero(p2.param)
            if is_phase_pi(p3) && p3.isparallel
                q3.param = zero(q3.param)
            else
                q3 = q3 + Phase(pi*im, true)
            end
        elseif is_phase_pi(p3) && p3.isparallel
            q1.param = zero(p1.param)
            q3.param = zero(p3.param)
            q2 = q2 + Phase(pi*im, true)
        end
    elseif is_phase_pi(p2) && p2.isparallel && is_phase_pi(p3) && p3.isparallel
        q2.param = zero(p1.param)
        q3.param = zero(p3.param)
        q1 = q1 + Phase(pi*im, true)
    end
    return q1, q2, q3
end
is_phase_pi(p::Phase) = (real(p.param) ≈ 0 && rem(imag(p.param)+pi, 2pi, RoundNearest) ≈ 0)

function yang_baxter_param(p1::Phase, p2::Phase, p3::Phase)
    # triangle => star
    p1.isparallel && (p1 = change_direction(p1))
    !p2.isparallel && (p2 = change_direction(p2))
    p3.isparallel && (p3 = change_direction(p3))

    q1_param, q2_param, q3_param = yang_baxter_param(p1.param, p2.param, p3.param)
    q1, q2, q3 = Phase(q1_param, true), Phase(q2_param, false), Phase(q3_param, true)
    return q1, q2, q3
end
function yang_baxter_param_inv(p1::Phase, p2::Phase, p3::Phase)
    # star => triangle
    !p1.isparallel && (p1 = change_direction(p1))
    p2.isparallel && (p2 = change_direction(p2))
    !p3.isparallel && (p3 = change_direction(p3))

    q1_param, q2_param, q3_param = yang_baxter_param_inv(p1.param, p2.param, p3.param)
    return Phase(q1_param, false), Phase(q2_param, true), Phase(q3_param, false)
end