using Test
using Quon
using Yao

function yb_fidelity(a1, b1, c1)
    circ1 = chain(1, Rz(-a1*im), Rx((-b1*im)), Rz(-c1*im))
    (a2, b2, c2) = yang_baxter_param(a1, b1, c1)
    circ2 = chain(1, Rx(-a2*im), Rz((-b2*im)), Rx(-c2*im))

    return operator_fidelity(matblock(mat(circ1') * mat(circ2)), I2)
end

a1 = rand()*π*im
b1 = rand()*π*im
c1 = rand()*π*im
a1, b1, c1 = rand(ComplexF64, 3)

a2, b2, c2 = yang_baxter_param(a1, b1, c1)
yb_fidelity(a1, b1, c1)

@test yb_fidelity(a1, b1, c1) ≈ 1
@test yb_fidelity(a2, b2, c2) ≈ 1

c1 = -a1
a2, b2, c2 = yang_baxter_param(a1, b1, c1)
@test yb_fidelity(a1, b1, c1) ≈ 1
@test yb_fidelity(a2, b2, c2) ≈ 1

c1 = π*im - a1
a2, b2, c2 = yang_baxter_param(a1, b1, c1)
@test yb_fidelity(a1, b1, c1) ≈ 1
@test yb_fidelity(a2, b2, c2) ≈ 1

@test_throws ErrorException yang_baxter_param(change_direction(pi/3*im), pi/3*im, change_direction(pi/3*im))
