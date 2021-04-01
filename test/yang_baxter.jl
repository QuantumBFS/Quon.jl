using Test

function check_inv(a, b, c)
    theta = [a,b,c] * im
    
    phi = yang_baxter_inv(yang_baxter(theta[1], theta[2], theta[3])...)
    phi = imag.([phi...])
    return (Rz(phi[3])*Rx(phi[2])*Rz(phi[1]))^-1 * Rz(c)*Rx(b)*Rz(a) ≈ [1 0; 0 1]
end
check_inv() = check_inv(rand(3)...)
Rz(θ) = exp(-θ/2*im*[1 0; 0 -1])
Rx(θ) = exp(-θ/2*im*[0 1; 1 0])

@test check_inv()