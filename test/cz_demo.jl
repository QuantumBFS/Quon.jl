using Quon
cz = tait_cz()
contract!(cz, tait_cz())
cz.locations[59] = (2.0, 1.0)
cz.locations[45] = (1.0, 1.0)
cz.locations[60] = (3.0, 1.0)
start_hes = Dict(
    1 => (3, -π/2),
    29 => (6, π),
    24 => (63, -π/2),
    45 => (26, 0),
    59 => (53, 0),
    60 => (54, π),
    57 => (131, π/4),
    32 => (87, π/2),
    56 => (153, π/2)
)
p = plot(cz; start_hes = start_hes, radius = 0.6, show_tait = true)
p |> SVG("1.svg")

m = match(Rule{:perm_rz}(), cz)
Quon.rewrite!(cz, m[3])
Quon.rewrite!(cz, m[9])
start_hes[45] = (4, π/2)
start_hes[60] = (64, π/2)
start_hes[57] = (103, π)
delete!(start_hes, 59)
p = plot(cz; radius = 0.5, start_hes = start_hes, show_tait = true)
p |> SVG("2.svg")

m_string_genus = match(Quon.Rule{:string_genus}(), cz)
Quon.rewrite!(cz, m_string_genus[1])
p = plot(cz; start_hes = start_hes, show_tait = true)
p |> SVG("3.svg")

m_fusion = match(Quon.Rule{:z_fusion}(), cz)
Quon.rewrite!(cz, m_fusion[1])
Quon.rewrite!(cz, m_fusion[2])
Quon.rewrite!(cz, m_fusion[3])
start_hes[57] = (156, -π/4)
p = plot(cz; start_hes = start_hes, show_tait = true)
p |> SVG("4.svg")

m_c_rm = match(Quon.Rule{:charge_rm_f}(), cz)
rewrite!(cz, m_c_rm[1])
start_hes[57] = (156, 0)
cz.locations[32] = (1.0, 3.0)
cz.locations[56] = (3.0, 3.0)
p = plot(cz; start_hes = start_hes, show_tait = true)
p |> SVG("5.svg")

m_gf = match(Rule(:genus_fusion), cz)
rewrite!(cz, m_gf[1])
start_hes[29] = (86, -3π/4)
cz.locations[29] = (2.0, 1.0)
p = plot(cz; start_hes = start_hes, show_tait = true)
p |> SVG("6.svg")
