using Quon, Test

cp = tait_copy()
contract_boundary_vertices!(cp, [1], [3])
@test plot(cp) !== nothing

rx0 = tait_rx(π*im)
contract_boundary_vertices!(rx0, [1], [6])
@test plot(Quon.simple_connected_planar_graph(rx0.g)[1]) !== nothing
@test plot(rx0) !== nothing

rz = tait_rz(im*pi)
contract!(rz, tait_rz(0.0*im))
@test plot(rz) !== nothing

contract_boundary_vertices!(rz, [1], [10])
@test plot(rz) !== nothing
rz = tait_rz(im*pi)
rx = tait_rx(im*pi)
contract!(rz, rx)
plot(rz)
@test plot(tensor_product!(rz, tait_id())) !== nothing

h = tait_hadamard()
id = tait_id()
c = tait_copy()
tensor_product!(id, h)
@test plot(id) !== nothing
contract!(c, id)
@test plot(c) !== nothing

c = contract!(tait_copy(), tait_copy(), [2, 3], [1, 3])
@test plot(c) !== nothing

s2 = contract!(tait_swap(), tait_swap())
s2.locations[26] = (1.0, 4.0)
s2.locations[29] = (3.0, 4.0)
s2.locations[25] = (2.0, 4.0)
s2.locations[33] = (2.0, 5.0)
s2.locations[30] = (1.0, 6.0)
s2.locations[34] = (2.0, 6.0)
s2.locations[32] = (3.0, 6.0)
s2.locations[27] = (1.0, 7.0)
s2.locations[31] = (2.0, 7.0)
s2.locations[28] = (3.0, 7.0)
s2.locations[19] = (0.0, 8.0)
s2.locations[23] = (2.0, 8.0)
s2.locations[20] = (4.0, 8.0)
st_hs = Dict(
    1 => (3, -π/4),
    8 => (50, -π/2),
    4 => (21, -3π/4),
    9 => (4, 3π/4),
    16 => (55, 0),
    12 => (22, π/4),
    13 => (29, 0),
    17 => (46, 0),
    15 => (41, 0),
    26 => (40, π/4),
    14 => (37, π/2),
    29 => (100, -π/4),
    25 => (34, π/2),
    22 => (26, π/4),
    24 => (24, π/4),
    33 => (111, 0),
    30 => (85, 0),
    34 => (102, 0),
    32 => (97, 0),
    27 => (66, -3π/4),
    28 => (72, -π/4),
    23 => (90, π/2),
    19 => (65, π/4),
    20 => (71, 3π/4)
)
@test plot(s2; backend = :compose, start_hes = st_hs) !== nothing
@test plot(s2; scale = 3) !== nothing