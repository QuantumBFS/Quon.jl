using Quon, Test

rx = tait_rx(π)
trace_vertex(rx, 1)
rz = tait_rz(π)
merge_graph!(rz, rx)
contract_boundary_vertices!(rz, [5], [6])

rz2 = tait_rz(im*pi)
rx2 = tait_rx(im*pi)
contract!(rz2, rx2)

@test nv(rz) == 6 == nv(rz2)
@test ne(rz) == 8 == ne(rz2)
@test nf(rz) == 4 == nf(rz2)
@test check_combinatorial_maps(rz) && check_combinatorial_maps(rz2)
@test check_vertices(rz) && check_faces(rz)
@test check_vertices(rz2) && check_faces(rz2)

idr = tensor_product!(tait_id(), rz2)
@test length(trace_vertex(idr, 4)) == 5

c = merge_graph!(tait_copy(), tait_copy())
contract_boundary_vertices!(c, [2], [9])
contract_boundary_vertices!(c, [3], [8])
@test c.g.vs_isolated[4] == 0

rx = tait_rx(im)
plot(rx)
rx.g.next[8]
rem_edge!(rx, 8; update = true)
rem_edge!(rx, 11; update = true)
@test length(rx.g.vs_isolated) == 1
