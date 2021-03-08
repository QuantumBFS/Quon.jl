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

