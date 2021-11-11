# TODO: handling singularity

const quon_atol = 1e-8

function change_direction(α::Complex)
    is_zero(α) && return InfZero
    is_pi(α) && return InfPi
    return to_quon_const(log(-tanh(α/2)))
end

"""
    yang_baxter_param(α1, β1, γ1)

Return α2, β2, γ2 according to the Yang-Baxter equation:

    \\    /    |       |    \\    / 
     \\  /     |       |     \\  /
      α1      |       |      α2
     |  \\    /         \\    /  |
     |   \\  /           \\  /   |
     |    β1      =      β2    |
     |   /  \\           /  \\   |
     |  /    \\         /    \\  |
      γ1      |       |      γ1 
     /  \\     |       |     /  \\
    /    \\    |       |    /    \\
"""
function yang_baxter_param(α1::Complex, β1::Complex, γ1::Complex)
    try
        if all(p -> (p isa QuonConst), to_quon_const.([α1, β1, γ1]))
            return yang_baxter_param(to_quon_const(α1), to_quon_const(β1), to_quon_const(γ1))
        end
        s1 = (α1+γ1)/2
        m1 = (α1-γ1)/2
        if !isapprox(cosh(s1), 0; atol = quon_atol)
            tanh_s2 = tanh(β1/2)*cosh(m1)/cosh(s1)
            @assert !(tanh_s2 ≈ 1)
            @assert !(tanh_s2 ≈ -1)
            s2 = atanh(tanh_s2)        
        elseif !isapprox(cosh(m1)*sinh(β1), 0; atol = quon_atol)
            s2 = π*im/2
        else
            error("This case should not be handled here")
        end
        if !isapprox(sinh(s1), 0; atol = quon_atol)
            tanh_m2 = -tanh(β1/2)*sinh(m1)/sinh(s1)
            @assert !(tanh_m2 ≈ 1)
            @assert !(tanh_m2 ≈ -1)
            m2 = atanh(tanh_m2)
        elseif !isapprox(sinh(m1)*sinh(β1), 0; atol = quon_atol)
            m2 = π*im/2
        else
            error("This case should not be handled here")
        end
        
        α2 = s2+m2
        γ2 = s2-m2
        if !isapprox(cosh(m2), 0; atol = quon_atol) && !isapprox(cosh(s1), 0; atol = quon_atol)
            β2 = 2*atanh(tanh(s1)*cosh(s2)/cosh(m2))
        elseif !isapprox(sinh(m2), 0; atol = quon_atol) && !isapprox(cosh(m1), 0; atol = quon_atol)
            β2 = 2*atanh(tanh(m1)*sinh(s2)/sinh(-m2))
        else
            error("This case should not be handled here")
        end
        
        return to_quon_const(α2), to_quon_const(β2), to_quon_const(γ2)
    catch e
        error("$((α1, β1, γ1)) is a singular point of Yang-Baxter equation")
    end
end

function yang_baxter_param_inv(α2, β2, γ2)
    a2 = change_direction(α2)
    b2 = change_direction(β2)
    α1, c1, b1 = yang_baxter_param(b2, a2, γ2)
    β1 = change_direction(b1)
    γ1 = change_direction(c1)
    return α1, β1, γ1
end
