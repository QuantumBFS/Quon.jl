using Quon

rz = tait_rz(im*pi)
contract_boundary_vertices!(rz, [1], [5])
plot(rz)
rz = tait_rz(im*pi)
rx = tait_rx(im*pi)
p1 = plot(rz; background = "white")
p2 = plot(rx; background = "white")
contract!(rz, rx)
plot(rz)
plot(tensor_product!(rz, tait_id()))

h = tait_hadamard()
id = tait_id()
c = tait_copy()
tensor_product!(id, h)
contract!(c, id)
plot(c)
tait_copy().outputs

c = contract!(tait_copy(), tait_copy(), [2, 3], [1, 3])

plot(c)

plot(tait_copy())

tait_copy().genuses
