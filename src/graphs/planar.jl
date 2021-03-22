# TODO: check the half edges are in reverse time order
struct PlanarGraph
    vertices::Dict{Vertex, HalfEdge}
    faces_half_edges::Dict{Face, HalfEdge} # face -> edge
    half_edges_faces::Dict{HalfEdge, Face} # half edge -> face
    nexts::Dict{HalfEdge, HalfEdge}
    prevs::Dict{HalfEdge, HalfEdge}
end

function PlanarGraph(half_edges_faces::Dict{HalfEdge, Face})
    faces_half_edges = Dict{Face, HalfEdge}()
    for (he, f) in half_edges_faces
        haskey(faces_half_edges, f) || (faces_half_edges[f] = he)
    end

    vertices = Dict{Vertex, HalfEdge}()
    for he in keys(half_edges_faces)
        haskey(vertices, he.src) || (vertices[he.src] = he)
    end
    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end

function PlanarGraph(vs::Dict{Vertex, HalfEdge}, f2he::Dict{Face, HalfEdge}, he2f::Dict{HalfEdge, Face})
    hes = keys(he2f)
    fs = Dict{Face, Vector{HalfEdge}}()
    for (he, f) in he2f
        if haskey(fs, f)
            push!(fs[f], he)
        else
            fs[f] = [he]
        end
    end
    nexts = Dict{HalfEdge, HalfEdge}()
    prevs = Dict{HalfEdge, HalfEdge}()
    for he in hes
        f_he = he2f[he]
        for nhe in fs[f_he]
            if nhe.src == he.dst
                nexts[he] = nhe
                prevs[nhe] = he
                break
            end
        end
    end
    return PlanarGraph(vs, f2he, he2f, nexts, prevs)
end

face(g::PlanarGraph, e::HalfEdge) = g.half_edges_faces[e]
nv(g::PlanarGraph) = length(keys(g.vertices))
src(g::PlanarGraph, e::HalfEdge) = src(e)
dst(g::PlanarGraph, e::HalfEdge) = dst(e)
twin(g::PlanarGraph, e::HalfEdge) = twin(e)

"""
    find_half_edge(g, v)

Find an out half edge of a vertex.
"""
find_half_edge(g::PlanarGraph, v::Vertex) = g.vertices[v]

"""
    find_half_edge(g, f)

Find a half edge around a face.
"""
find_half_edge(g::PlanarGraph, f::Face) = g.faces_half_edges[f]

"""
    next(g, he)

Returns the next half edge of a half edge.
"""
next(g::PlanarGraph, e::HalfEdge) = g.nexts[e]

"""
    prev(g, he)

Returns the previous half edge of a half edge.
"""
prev(g::PlanarGraph, e::HalfEdge) = g.prevs[e]
σ(g::PlanarGraph, e::HalfEdge) = twin(g, prev(g, e))
σ_inv(g::PlanarGraph, e::HalfEdge) = next(g, twin(g, e))

has_half_edge(g::PlanarGraph, he::HalfEdge) = haskey(g.half_edges_faces, he)
has_half_edge(g::PlanarGraph, s::Integer, d::Integer) = has_half_edge(g, HalfEdge(s, d))
has_vertex(g::PlanarGraph, v::Vertex) = haskey(g.vertices, v)
has_vertex(g::PlanarGraph, v::Integer) = has_vertex(g, Vertex(v))

"""
    is_boundary(g, he)

Check whether a half edge `he` is on the boundary of a planar graph `g`.
"""
is_boundary(g::PlanarGraph, he::HalfEdge) = face(g, he) == Face(0)

"""
    trace_face(g, f)

Returns a tuple containing all half edges and all vertices around a face.
"""
function trace_face(g::PlanarGraph, f::Face)
    starting_he = find_half_edge(g, f)
    hes = [starting_he]
    vs = [src(starting_he)]

    he = next(g, starting_he)
    while he != starting_he
        push!(hes, he)
        push!(vs, src(he))
        he = next(g, he)
    end

    return hes, vs
end
trace_face(g::PlanarGraph, f::Integer) = trace_face(g, Face(f))

"""
    trace_vertex(g, v)

Returns a tuple containing all out half edges and all out neighbors from a vertex.
"""
function trace_vertex(g::PlanarGraph, v::Vertex)
    starting_he = find_half_edge(g, v)
    hes = [starting_he]
    out_nbs = [dst(starting_he)]

    he = σ(g, starting_he)
    while he != starting_he
        push!(hes, he)
        push!(out_nbs, dst(he))
        he = σ(g, he)
    end

    return hes, out_nbs
end
trace_vertex(g::PlanarGraph, v::Integer) = trace_vertex(g, Vertex(v))

"""
    vertices(g)

Returns all vertices in a planar graph.
"""
vertices(g::PlanarGraph) = collect(keys(g.vertices))

faces(g::PlanarGraph) = collect(keys(g.faces_half_edges))

"""
    half_edges(g)

Returns all half edges in a planar graph `g`.
"""
half_edges(g::PlanarGraph) = collect(keys(g.half_edges_faces))

"""
    edges(g)

Return all whole edges in a planar graph `g`.
"""
function edges(g::PlanarGraph)
    es = Tuple{Vertex, Vertex}[]
    hes = half_edges(g)
    for he in hes
        if src(he).id < dst(he).id
            push!(es, (he.src, he.dst))
        end
    end
    return es
end

"""
    planar_rz()

Generate the following planar graph for Rz:

        1 ------- 2 ----- 3 ----- 4
       / \\\\ (1)  //      ||     ||  \\
      /   \\\\    //       ||     ||   \\
     /     \\\\  //  (2)   || (3) ||    \\
    9  (6)   10          11     12 (4) 13
     \\     //  \\\\        ||     ||    /
      \\   //    \\\\       ||     ||   / 
       \\ // (5)  \\\\      ||     ||  /
        8 ------- 7 ----- 6 ----- 5
"""
function planar_rz()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 10),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 13),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 8),
        Vertex(10) => HalfEdge(10, 2),
        Vertex(11) => HalfEdge(11, 3),
        Vertex(12) => HalfEdge(12, 4),
        Vertex(13) => HalfEdge(13, 4),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(0) => HalfEdge(1, 2),
        Face(1) => HalfEdge(1, 10),
        Face(2) => HalfEdge(2, 10),
        Face(3) => HalfEdge(3, 11),
        Face(4) => HalfEdge(4, 12),
        Face(5) => HalfEdge(8, 7),
        Face(6) => HalfEdge(1, 9),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(1, 2) => Face(0),
        HalfEdge(2, 3) => Face(0),
        HalfEdge(3, 4) => Face(0),
        HalfEdge(4, 13) => Face(0),
        HalfEdge(13, 5) => Face(0),
        HalfEdge(5, 6) => Face(0),
        HalfEdge(6, 7) => Face(0),
        HalfEdge(7, 8) => Face(0),
        HalfEdge(8, 9) => Face(0),
        HalfEdge(9, 1) => Face(0),

        HalfEdge(2, 1) => Face(1),
        HalfEdge(1, 10) => Face(1),
        HalfEdge(10, 2) => Face(1),

        HalfEdge(2, 10) => Face(2),
        HalfEdge(10, 7) => Face(2),
        HalfEdge(7, 6) => Face(2),
        HalfEdge(6, 11) => Face(2),
        HalfEdge(11, 3) => Face(2),
        HalfEdge(3, 2) => Face(2),

        HalfEdge(3, 11) => Face(3),
        HalfEdge(11, 6) => Face(3),
        HalfEdge(6, 5) => Face(3),
        HalfEdge(5, 12) => Face(3),
        HalfEdge(12, 4) => Face(3),
        HalfEdge(4, 3) => Face(3),

        HalfEdge(4, 12) => Face(4),
        HalfEdge(12, 5) => Face(4),
        HalfEdge(5, 13) => Face(4),
        HalfEdge(13, 4) => Face(4),

        HalfEdge(8, 7) => Face(5),
        HalfEdge(7, 10) => Face(5),
        HalfEdge(10, 8) => Face(5),

        HalfEdge(1, 9) => Face(6),
        HalfEdge(9, 8) => Face(6),
        HalfEdge(8, 10) => Face(6),
        HalfEdge(10, 1) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end

"""
    planar_rx()

Generate the following planar graph for Rx:

         1 ----- 2 ------ 3 ----- 4
       /  ||     \\\\ (3)  //     ||  \\
      /   ||      \\\\    //      ||   \\
     /    ||  (2)  \\\\  //  (4)  ||    \\
    9 (1) 10         11         12 (5) 13
     \\    ||       //  \\\\       ||    /
      \\   ||      //    \\\\      ||   / 
       \\  ||     // (6)  \\\\     ||  /
         8 ----- 7 ------ 6 ----- 5
"""
function planar_rx()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 10),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 13),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 8),
        Vertex(10) => HalfEdge(10, 1),
        Vertex(11) => HalfEdge(11, 2),
        Vertex(12) => HalfEdge(12, 4),
        Vertex(13) => HalfEdge(13, 4),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(0) => HalfEdge(1, 2),
        Face(1) => HalfEdge(1, 9),
        Face(2) => HalfEdge(2, 1),
        Face(3) => HalfEdge(3, 2),
        Face(4) => HalfEdge(4, 3),
        Face(5) => HalfEdge(5, 13),
        Face(6) => HalfEdge(6, 11),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(1, 2) => Face(0),
        HalfEdge(2, 3) => Face(0),
        HalfEdge(3, 4) => Face(0),
        HalfEdge(4, 13) => Face(0),
        HalfEdge(13, 5) => Face(0),
        HalfEdge(5, 6) => Face(0),
        HalfEdge(6, 7) => Face(0),
        HalfEdge(7, 8) => Face(0),
        HalfEdge(8, 9) => Face(0),
        HalfEdge(9, 1) => Face(0),

        HalfEdge(1, 9) => Face(1),
        HalfEdge(9, 8) => Face(1),
        HalfEdge(8, 10) => Face(1),
        HalfEdge(10, 1) => Face(1),

        HalfEdge(2, 1) => Face(2),
        HalfEdge(1, 10) => Face(2),
        HalfEdge(10, 8) => Face(2),
        HalfEdge(8, 7) => Face(2),
        HalfEdge(7, 11) => Face(2),
        HalfEdge(11, 2) => Face(2),

        HalfEdge(2, 11) => Face(3),
        HalfEdge(11, 3) => Face(3),
        HalfEdge(3, 2) => Face(3),

        HalfEdge(4, 3) => Face(4),
        HalfEdge(3, 11) => Face(4),
        HalfEdge(11, 6) => Face(4),
        HalfEdge(6, 5) => Face(4),
        HalfEdge(5, 12) => Face(4),
        HalfEdge(12, 4) => Face(4),

        HalfEdge(5, 13) => Face(5),
        HalfEdge(13, 4) => Face(5),
        HalfEdge(4, 12) => Face(5),
        HalfEdge(12, 5) => Face(5),

        HalfEdge(6, 11) => Face(6),
        HalfEdge(11, 7) => Face(6),
        HalfEdge(7, 6) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end
