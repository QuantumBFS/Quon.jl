using Test
using Quon
using Yao

function check_inv(a, b, c)
    theta = [a,b,c] * im
    
    phi = yang_baxter_param_inv(yang_baxter_param(theta[1], theta[2], theta[3])...)
    phi = imag.([phi...])

    circ1 = chain(1, Rz(a), Rx(b), Rz(c))
    circ2 = chain(1, Rz(phi[1]), Rx(phi[2]), Rz(phi[3]))
    return operator_fidelity(circ1, circ2) ≈ 1
end
check_inv() = check_inv(rand(3)...)
check_inv()

@test check_inv()