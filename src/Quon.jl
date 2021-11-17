module Quon

export yang_baxter_param, yang_baxter_param_inv, change_direction, change_direction!
export PlanarMultigraph, Tait,
    vertices, faces, half_edges,
    src, dst, half_edge, face,
    next, prev, twin, α, ϕ, σ, σ_inv,
    nv, ne, nf, nhe, is_boundary,
    out_half_edge, surrounding_half_edge,
    trace_face, trace_vertex, neighbors,
    rem_vertex!, rem_edge!, rem_face!, merge_graph!,
    add_edge!, add_vertex!,
    check_faces, check_vertices, check_combinatorial_maps
export contract!, contract_boundary_vertices!, 
    tensor_product!, merge_boundary_vertices!,
    quon_param, quon_params, genuses, is_genus, is_open_vertex
export planar_rx, planar_rz, planar_id, planar_copy, planar_zero_state, planar_swap,
    tait_rx, tait_rz, tait_id, tait_hadamard, tait_copy, tait_cz, tait_cnot, 
    tait_zero_state, tait_swap

include("utils/yang_baxter.jl")

include("params/quon_const.jl")
include("params/quon_param.jl")

include("graphs/planar_multigraph.jl")
include("graphs/tait.jl")
include("graphs/predefined.jl")

export Rule, rewrite!, match!

include("match.jl")
include("rewrite.jl")
include("check.jl")
include("simplify.jl")

include("utils/zx.jl")
include("circuits.jl")

export push_gate!

export plot

include("plots.jl")

end
