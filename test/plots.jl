using Quon

rz = tait_rz(im*pi)
rx = tait_rx(im*pi)
p1 = plot(rz; show_faces = false, show_half_edges = false, background = "white")
# p1 |> SVG("rz.svg")
p2 = plot(rx; background = "white")
# p2 |> SVG("rzrx_all.svg")
contract!(rz, rx)
plot(rz)
id = tait_id()
plot(contract!(tait_id(), tait_id()))
