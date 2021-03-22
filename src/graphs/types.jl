# 1. charge is encoded in the vertices of the medial graph
# 2. 

# import Base: ==

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
src(e::HalfEdge) = e.src
dst(e::HalfEdge) = e.dst
twin(e::HalfEdge) = HalfEdge(e.dst, e.src)

# ==(he1::HalfEdge, he2::HalfEdge) = (he1.src == he2.src && he1.dst == he2.dst)