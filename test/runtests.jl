using Quon
using Test

Quon.planar_rz()

vertices = Dict{Vertex, HalfEdge}(
    Vertex(1) => HalfEdge(1, 5),
    Vertex(2) => HalfEdge(2, 5),
    Vertex(3) => HalfEdge(3, 5),
    Vertex(4) => HalfEdge(4, 5),
    Vertex(5) => HalfEdge(5, 3),
)

faces = Dict{Face, HalfEdge}(
    Face(1) => HalfEdge(1, 5),
    Face(2) => HalfEdge(2, 5),
    Face(3) => HalfEdge(4, 5),
    Face(4) => HalfEdge(3, 5),
)

halfedge_faces = Dict{HalfEdge, Face}(
    HalfEdge(1, 5) => Face(1),
    HalfEdge(5, 2) => Face(1),
    HalfEdge(2, 5) => Face(2),
    HalfEdge(5, 4) => Face(2),
    HalfEdge(4, 5) => Face(3),
    HalfEdge(5, 3) => Face(3),
    HalfEdge(3, 5) => Face(4),
    HalfEdge(5, 1) => Face(4),
)

nexts = Dict{HalfEdge, HalfEdge}(
    HalfEdge(1, 5) => HalfEdge(5, 2),
    HalfEdge(2, 5) => HalfEdge(5, 4),
    HalfEdge(4, 5) => HalfEdge(5, 3),
    HalfEdge(3, 5) => HalfEdge(5, 1),
)

prevs = Dict{HalfEdge, HalfEdge}(
    HalfEdge(5, 2) => HalfEdge(1, 5),
    HalfEdge(5, 4) => HalfEdge(2, 5),
    HalfEdge(5, 3) => HalfEdge(4, 5),
    HalfEdge(5, 1) => HalfEdge(3, 5),
)

graph = PlanarGraph(vertices, faces, HalfEdgeTable(halfedge_faces, nexts, prevs))
