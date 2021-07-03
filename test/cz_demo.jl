using Quon, Compose
cz = tait_cz()
contract!(cz, tait_cz())
cz.locations[4] = (2.0, 1.0)
start_hes = Dict(
    1 => (3, -π/2),
    5 => (6, π),
    24 => (63, -π/2),
    7 => (26, 0),
    4 => (53, 0),
    16 => (54, π),
    34 => (131, π/4),
    32 => (87, π/2),
    56 => (153, π/2)
)
p = plot(cz; start_hes = start_hes, radius = 0.6, show_tait = false)
p |> SVG("1.svg")

m = match(Rule{:perm_rz}(), cz)
Quon.rewrite!(cz, m[3])
Quon.rewrite!(cz, m[9])
start_hes[7] = (4, π/2)
start_hes[16] = (64, π/2)
start_hes[34] = (103, π)
delete!(start_hes, 4)
p = plot(cz; radius = 0.5, start_hes = start_hes, show_tait = false)
p |> SVG("2.svg")

m_string_genus = match(Quon.Rule{:string_genus}(), cz)
Quon.rewrite!(cz, m_string_genus[1])
p = plot(cz; start_hes = start_hes, show_tait = false)
p |> SVG("3.svg")

m_fusion = match(Quon.Rule{:z_fusion}(), cz)
Quon.rewrite!(cz, m_fusion[1])
Quon.rewrite!(cz, m_fusion[2])
Quon.rewrite!(cz, m_fusion[3])
start_hes[34] = (158, π/4)
p = plot(cz; start_hes = start_hes, show_tait = false)
p |> SVG("4.svg")

m_yb = match(Quon.Rule{:yang_baxter_triangle}(), cz)
rewrite!(cz, m_yb[1])
p = plot(cz; start_hes = start_hes, show_tait = false)
p |> SVG("5.svg")

m_id = match(Quon.Rule{:identity}(), cz)
Quon.rewrite!(cz, m_id[1])
Quon.rewrite!(cz, m_id[3])
Quon.rewrite!(cz, m_id[5])
start_hes[34] = (86, π)
cz.locations[32] = (1.0, 3.0)
cz.locations[56] = (3.0, 3.0)
p = plot(cz; start_hes = start_hes, show_tait = false)
p |> SVG("6.svg")

m_gf = match(Rule(:genus_fusion), cz)
rewrite!(cz, m_gf[2])
start_hes[34] = (86, -3π/4)
cz.locations[34] = (2.0, 1.0)
p = plot(cz; start_hes = start_hes, show_tait = false)
p |> SVG("7.svg")

using Compose
compose(context(),
	bezigon((0, 0), [[(0, 1), (0, 1), (1, 1)], [(1, 0), (1, 0), (0, 0)]])
)