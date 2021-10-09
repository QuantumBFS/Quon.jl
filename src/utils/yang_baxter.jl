# TODO: handling singularity

const quon_atol = 1e-8

change_direction(α) = log(-tanh(α/2))

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
function yang_baxter_param(α1, β1, γ1)
    s1 = (α1+γ1)/2
    m1 = (α1-γ1)/2
    if !isapprox(cosh(s1), 0; atol = quon_atol)
        tanh_s2 = tanh(β1/2)*cosh(m1)/cosh(s1)
        @assert !(tanh_s2 ≈ 1) && !(tanh_s2 ≈ -1)
        s2 = atanh(tanh_s2)
    elseif !isapprox(cosh(m1), 0; atol = quon_atol)
        s2 = π*im/2
    else
        s2 = β1/2
    end
    if !isapprox(sinh(s1), 0; atol = quon_atol)
        tanh_m2 = -tanh(β1/2)*sinh(m1)/sinh(s1)
        @assert !(tanh_m2 ≈ 1) && !(tanh_m2 ≈ -1)
        m2 = atanh(tanh_m2)
    elseif !isapprox(sinh(m1), 0; atol = quon_atol)
        m2 = π*im/2
    else
        m2 = -β1/2
    end
    
    α2 = s2+m2
    γ2 = s2-m2
    if !isapprox(cosh(m2), 0; atol = quon_atol)
        β2 = 2*atanh(tanh(s1)*cosh(s2)/cosh(m2))
    elseif isapprox(sinh(s1), 0; atol = quon_atol)
        tanh_β2_2 = -cosh(s2)/cosh(s1)*sinh(m1)/sinh(m2)*tanh(β1/2)
        β2 = 2*atanh(tanh_β2_2)
    elseif !isapprox(cosh(s2), 0; atol = quon_atol)
        β2 = pi*im
    else
        β2 = 2*s1
    end
    
    return α2, β2, γ2
end

function yang_baxter_param_inv(α2, β2, γ2)
    a2 = change_direction(α2)
    b2 = change_direction(β2)
    α1, c1, b1 = yang_baxter_param(b2, a2, γ2)
    β1 = change_direction(b1)
    γ1 = change_direction(c1)
    return α1, β1, γ1
end
