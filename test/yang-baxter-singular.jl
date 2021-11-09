using Test
using Quon
using Yao

function yb_fidelity(a1, b1, c1)
    circ1 = chain(1, Rz(-a1*im), Rx((-b1*im)), Rz(-c1*im))
    (a2, b2, c2) = yang_baxter_param(a1, b1, c1)
    circ2 = chain(1, Rx(-a2*im), Rz((-b2*im)), Rx(-c2*im))

    return operator_fidelity(circ1, circ2)
end

a = rand()*π*im
b = rand()*π*im
c = rand()*π*im

@test yb_fidelity(a, b, c) ≈ 1
@test yb_fidelity(yang_baxter_param(a, b, c)...) ≈ 1

c = -a
@test yb_fidelity(a, b, c) ≈ 1
@test yb_fidelity(yang_baxter_param(a, b, c)...) ≈ 1

c = π*im - a
@test yb_fidelity(a, b, c) ≈ 1
@test yb_fidelity(yang_baxter_param(a, b, c)...) ≈ 1

@test_throws ErrorException yang_baxter_param(pi/2*im, pi/2*im, pi/2*im)
@test_throws ErrorException yang_baxter_param(change_direction(pi/3*im), pi/3*im, change_direction(pi/3*im))
