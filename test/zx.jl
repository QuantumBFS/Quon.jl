using Quon
using ZXCalculus

q = Tait(2)
push_gate!(q, Val(:CZ), 1, 2)
push_gate!(q, Val(:CZ), 1, 2)
plot(q)
plot(q; backend = :compose)
zxg = ZXDiagram(q) |> ZXGraph |> full_reduction
zxg |> ancilla_extraction |> plot

using ZXCalculus
# s = tait_swap()
# zxd = ZXDiagram(s)
# zxg = ZXGraph(zxd)
# full_reduction(zxg)
# circ = ancilla_extraction(zxg) 
# circ |> plot

# zxg = ZXGraph(zxd)
# YaoPlots.plot(zxg)
# clifford_simplification(zxg)
# YaoPlots.plot(zxg)
# ZXCalculus.generate_layout!(zxg)

# rz = tait_rz(pi*im)
# zxd_rz = ZXDiagram(rz)
# zxd_h = ZXDiagram(tait_hadamard())
# zxd_cz = ZXDiagram(tait_cz())
zxg_cnot(i, j) = tait_cnot(i, j) |> ZXDiagram |> ZXGraph |> full_reduction |> ancilla_extraction
zxg_cnot(1, 30) |> plot