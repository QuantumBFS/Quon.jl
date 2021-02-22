# 1. charge is encoded in the vertices of the medial graph
# 2. 

struct Vertex
    id::Int
end

struct Face
    id::Int
end

Base.convert(::Type{Vertex}, id::Int) = Vertex(id)
Base.convert(::Type{Face}, id::Int) = Face(id)

struct HalfEdge
    src::Vertex
    dst::Vertex
end

HalfEdge(src::Int, dst::Int) = HalfEdge(Vertex(src), Vertex(dst))

struct PlanarGraph
    vertices::Dict{Vertex, HalfEdge}
    faces_half_edges::Dict{Face, HalfEdge} # face -> edge
    half_edges_faces::Dict{HalfEdge, Face} # half edge -> face
    # nexts::Dict{HalfEdge, HalfEdge}
    # prevs::Dict{HalfEdge, HalfEdge}
end

face(g::PlanarGraph, e::HalfEdge) = g.half_edges_faces[e]

function planar_rz()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 9),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 4),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 2),
        Vertex(10) => HalfEdge(10, 3),
        Vertex(11) => HalfEdge(11, 4),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(1) => HalfEdge(1, 9),
        Face(2) => HalfEdge(2, 9),
        Face(3) => HalfEdge(3, 10),
        Face(4) => HalfEdge(4, 11),
        Face(5) => HalfEdge(8, 7),
        Face(6) => HalfEdge(1, 8),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(2, 1) => Face(1),
        HalfEdge(1, 9) => Face(1),
        HalfEdge(9, 2) => Face(1),

        HalfEdge(2, 9) => Face(2),
        HalfEdge(9, 7) => Face(2),
        HalfEdge(7, 6) => Face(2),
        HalfEdge(6, 10) => Face(2),
        HalfEdge(10, 3) => Face(2),
        HalfEdge(3, 2) => Face(2),

        HalfEdge(3, 10) => Face(3),
        HalfEdge(10, 6) => Face(3),
        HalfEdge(6, 5) => Face(3),
        HalfEdge(5, 11) => Face(3),
        HalfEdge(11, 4) => Face(3),
        HalfEdge(4, 3) => Face(3),

        HalfEdge(4, 11) => Face(4),
        HalfEdge(11, 5) => Face(4),
        HalfEdge(5, 4) => Face(4),

        HalfEdge(8, 7) => Face(5),
        HalfEdge(7, 9) => Face(5),
        HalfEdge(9, 8) => Face(5),

        HalfEdge(1, 8) => Face(6),
        HalfEdge(8, 9) => Face(6),
        HalfEdge(9, 1) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end

function planar_rx()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 9),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 4),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 8),
        Vertex(10) => HalfEdge(10, 2),
        Vertex(11) => HalfEdge(11, 5),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(1) => HalfEdge(1, 8),
        Face(2) => HalfEdge(2, 1),
        Face(3) => HalfEdge(3, 2),
        Face(4) => HalfEdge(7, 6),
        Face(5) => HalfEdge(4, 3),
        Face(6) => HalfEdge(5, 4),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(9, 1) => Face(1),
        HalfEdge(1, 8) => Face(1),
        HalfEdge(8, 9) => Face(1),

        HalfEdge(1, 9) => Face(2),
        HalfEdge(9, 8) => Face(2),
        HalfEdge(8, 7) => Face(2),
        HalfEdge(7, 10) => Face(2),
        HalfEdge(10, 2) => Face(2),
        HalfEdge(2, 1) => Face(2),

        HalfEdge(2, 10) => Face(3),
        HalfEdge(10, 3) => Face(3),
        HalfEdge(3, 2) => Face(3),

        HalfEdge(6, 10) => Face(4),
        HalfEdge(10, 7) => Face(4),
        HalfEdge(7, 6) => Face(4),

        HalfEdge(3, 10) => Face(5),
        HalfEdge(10, 6) => Face(5),
        HalfEdge(6, 5) => Face(5),
        HalfEdge(5, 11) => Face(5),
        HalfEdge(11, 4) => Face(5),
        HalfEdge(4, 3) => Face(5),

        HalfEdge(4, 11) => Face(6),
        HalfEdge(11, 5) => Face(6),
        HalfEdge(5, 4) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end
