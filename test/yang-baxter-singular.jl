using Quon

yang_baxter_param(pi/2*im, pi*im, pi/2*im)

using Quon: QuonParam
p1, p2, p3 = QuonParam(rand()*2π*im, true), QuonParam(rand()*2π*im, false), QuonParam(rand()*2π*im, true)
q1, q2, q3 = yang_baxter_param(p1, p2, p3)
change_direction.([q1, q2, q3])
(collect(yang_baxter_param(q1, q2, q3)))
cdr = change_direction
a, b = (cdr(q1), cdr(q2))
r1, r2, r3 = yang_baxter_param(QuonParam(π/4*im, false), QuonParam(π/4*im, true), QuonParam(π/4*im, false))
cdr(r1).param

a, b, c = yang_baxter_param(3π/4*im, π/2*im, 5π/4*im)
(d, e, f) = yang_baxter_param(a, b, c)
p, q, r = yang_baxter_param(π/4*im, b, π/4*im)
c+r
yang_baxter_param(a+p, q, c+r)

t = [3π/4*im, π/2*im, 5π/4*im]
yb_eq(v) = [yang_baxter_param(v[1], v[2], v[3])...]
yb_eq2(v) = yb_eq(yb_eq(v))

isapprox(yb_eq2(t), t)


Rz(x) = exp(-x/2*[1 0; 0 -1])
Rx(x) = exp(-x/2*[0 1; 1 0])

p1, p2, p3 = rand(3)*2π*im
q1, q2, q3 = yang_baxter_param(p1, p2, p3)
r2, r3, r1 = yang_baxter_param(p3, change_direction(p1), change_direction(p2))
r2 = change_direction(r2)
r3 = change_direction(r3)

MP = Rz(p3) * Rx(p2) * Rz(p1)
MQ = Rx(q3) * Rz(q2) * Rx(q1)
MR = Rx(r3) * Rz(r2) * Rx(r1)
MR * MQ'

M1 = Rz(5π/4)*Rx(π/2)*Rz(3π/4)
M2 = Rx(-atan(sqrt(2)/2)-π/2)*
    Rz(2*atan(-1/sqrt(3)))*
    Rx(-atan(sqrt(2)/2)+π/2)
M1' * M2

t = rand(3)*2pi*im
collect(yang_baxter_param_inv(t...)) - collect(yang_baxter_param(t...))
yb_eq2(t)

yb_eq(im*[π/2, 5π/4, -π/2])