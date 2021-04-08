using Quon, Test

rx = tait_rx(π)
trace_vertex(rx, 1)
rz = tait_rz(π)
merge_graph!(rz, rx)
contract_boundary_vertices!(rz, [5], [6])
@test nv(rz) == 6
@test ne(rz) == 8
@test nf(rz) == 4
@test check_combinatorial_maps(rz)
