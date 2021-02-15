# 1. charge is encoded in the vertices of the medial graph
# 2. 

struct Vertex
    id::Int
    coordinate::NTuple{2, Int}
end

struct Face
    id::Int
end

struct HalfEdge
    id::Int
    src::Vertex
    dst::Vertex
end

struct IncidentVertex
    v::Vertex
    incident_edge::Int
end

struct HalfEdgeTable
    twins::Dict{Int, HalfEdge} # half edge -> twin id
    faces::Dict{Int, Face} # half edge -> face
    nexts::Dict{Int, HalfEdge}
    prevs::Dict{Int, HalfEdge}
end

struct PlanarGraph
    vertices::Dict{Int, Vertex}
    faces::Dict{Face, HalfEdge} # face -> edge
    he_table::HalfEdgeTable
end

twin(g::PlanarGraph, e::HalfEdge) = g.he_table.twins[e.id]
face(g::PlanarGraph, e::HalfEdge) = g.he_table.faces[e.id]
next_half_edge(g::PlanarGraph, e::HalfEdge) = g.he_table.nexts[e.id]

function half_edges(g::PlanarGraph, f::Face)
    start = g.faces[f]
    edges = HalfEdge[start]
    curr = next_half_edge(g, start)

    while curr != start
        push!(edges, curr)
        curr = next_half_edge(g, curr)
    end
    return edges
end
