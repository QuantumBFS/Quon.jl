using Quon
circ = tait_swap()
contract!(circ, tait_swap())
# plot(circ)

Quon.simplify!(circ, Quon.Rule(:swap_genus))
match(Quon.Rule(:x_fusion), circ)
Quon.simplify!(circ, Quon.Rule(:x_fusion))
match(Quon.Rule(:yang_baxter_triangle), circ)

plot(circ; backend = :compose)

using ZXCalculus
zxg = ZXGraph(ZXDiagram(circ))
zxg |> full_reduction |> ancilla_extraction |> plot