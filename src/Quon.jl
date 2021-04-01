module Quon

include("utils/yang_baxter.jl")
export yang_baxter, yang_baxter_inv, change_direction

# export PlanarGraph, Vertex, Face, HalfEdge, 
#     # interfaces
#     half_edges, find_half_edge,
#     twin, nv, src, dst, twin, next, prev, face, σ, σ_inv,
#     has_vertex, has_half_edge, 
#     is_boundary, trace_face, trace_vertex, 
#     contract,
#     half_edges, vertices,
#     edges,

#     # predefined
#     planar_rx,
#     planar_rz

# export QuonGraph, quon_rx, quon_rz

include("graphs/tait.jl")
export TaitGraph, Edge,
    # APIs
    nv, ne, src, dst, edge, neighbors,
    rem_edge!, rem_vertex!,
    another_end, adjacent_edges, 
    contract!, contract_boundary_vertices!,
    
    # predefined
    tait_rx, tait_rz 

# include("graphs/types.jl")
# include("graphs/planar.jl")
# include("quon_graph.jl")
# include("contract.jl")
# include("plots.jl")

end
