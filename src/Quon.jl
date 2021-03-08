module Quon

export PlanarGraph, HalfEdgeTable, Vertex, Face, HalfEdge,
    # interfaces
    half_edges,
    twin,
    face,
    contract,
    vertices_map,
    faces_map,
    # predefined
    planar_rx,
    planar_rz

include("graphs/types.jl")
include("graphs/planar.jl")
include("graphs/medial.jl")
include("contract.jl")

end
