module Quon

export yang_baxter_param, yang_baxter_param_inv, change_direction!
export PlanarMultigraph, QuonTait,
    vertices, faces, half_edges,
    src, dst, half_edge, face,
    next, prev, twin, α, ϕ, σ, σ_inv,
    nv, ne, nf, nhe, is_boundary,
    out_half_edge, surrounding_half_edge,
    trace_face, trace_vertex, neighbors,
    rem_vertex!, rem_edge!, rem_face!, merge_graph!,
    check_faces, check_vertices, check_combinatorial_maps
export contract!, contract_boundary_vertices!, 
    tensor_product!, merge_boundary_vertices!,
    phases, phases, genuses, is_genus, is_open
export plot
export planar_rx, planar_rz, planar_id, tait_rx, tait_rz, tait_id 

include("utils/yang_baxter.jl")
include("graphs/planar_multigraph.jl")
include("graphs/tait.jl")
include("graphs/predefined.jl")
include("plots.jl")

end
