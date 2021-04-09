module Quon

include("utils/yang_baxter.jl")
export yang_baxter_param, yang_baxter_param_inv, change_direction!

include("graphs/planar_multigraph.jl")
export PlanarMultigraph, QuonTait,
    vertices, faces, half_edges,
    src, dst, half_edge, face,
    next, prev, twin, α, ϕ, σ, σ_inv,
    nv, ne, nf, nhe, is_boundary,
    out_half_edge, surrounding_half_edge,
    trace_face, trace_vertex, neighbors,
    rem_vertex!, rem_edge!, rem_face!, merge_graph!,
    check_faces, check_vertices, check_combinatorial_maps

include("graphs/tait.jl")
export contract!, contract_boundary_vertices!, 
    phases, phases, genuses, is_genus, is_open

include("graphs/predefined.jl")
export planar_rx, planar_rz, tait_rx, tait_rz 

include("plots.jl")
export plot

end
