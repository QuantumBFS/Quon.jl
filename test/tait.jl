using Quon
using Test

rx = tait_rx(π)
rz = tait_rz(π)

contract!(rz, rx, [5], [5])
@test nv(rz) == 6
@test ne(rz) == 8
@test rz.boundaries == [4, 11]