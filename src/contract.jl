function check_indices(a::PlanarGraph, b::PlanarGraph, Ia::Vector{Vertex}, Ib::Vector{Vertex})
    @assert length(Ia) == length(Ib) "length of indices does not match"
    if !haskey(a.half_edges_faces, HalfEdge(Ia[1], Ia[2]))
        Ia = reverse(Ia)
        Ib = reverse(Ib)
    end

    n = length(Ia)
    for k in 1:n-1
        he_a = HalfEdge(Ia[k], Ia[k+1])
        he_b = HalfEdge(Ib[k+1], Ib[k])
        he_a in keys(a.half_edges_faces) && he_b in keys(b.half_edges_faces) ||
            error("invalid contraction indices, the indices" *
                " order might be incorrect, or contains non-open edge")
    end
    return Ia, Ib
end

function vertices_map(a::PlanarGraph, b::PlanarGraph, Ia::Vector{Vertex}, Ib::Vector{Vertex})
    vertices_map_a = Dict{Vertex, Vertex}() # original -> new
    vertices_map_b = Dict{Vertex, Vertex}() # original -> new

    n_vertices = 1
    for k in 1:length(Ia)
        v_a = Ia[k]
        v_b = Ib[k]
        vertices_map_a[v_a] = Vertex(n_vertices)
        vertices_map_b[v_b] = Vertex(n_vertices)
        n_vertices += 1
    end

    for v in sort!(collect(keys(a.vertices)), by = x -> x.id)
        if !haskey(vertices_map_a, v)
            vertices_map_a[v] = Vertex(n_vertices)
            n_vertices += 1
        end
    end

    for v in sort!(collect(keys(b.vertices)), by = x -> x.id)
        if !haskey(vertices_map_b, v)
            vertices_map_b[v] = Vertex(n_vertices)
            n_vertices += 1
        end
    end

    return vertices_map_a, vertices_map_b
end

function faces_map(a::PlanarGraph, b::PlanarGraph, Ia::Vector{Vertex}, Ib::Vector{Vertex})
    faces_map_a = Dict{Face, Face}()
    faces_map_b = Dict{Face, Face}() 
    
    n_face = 1
    for k in 1:length(Ia)-1
        he_a = HalfEdge(Ia[k], Ia[k+1])
        he_b = HalfEdge(Ib[k+1], Ib[k])
        fa = a.half_edges_faces[he_a]
        fb = b.half_edges_faces[he_b]
        # give f a new id
        faces_map_a[fa] = Face(n_face)
        faces_map_b[fb] = Face(n_face)
        n_face += 1
    end

    for f in sort!(collect(keys(a.faces_half_edges)), by = x -> x.id)
        if !haskey(faces_map_a, f)
            faces_map_a[f] = Face(n_face)
            n_face += 1
        end
    end

    for f in sort!(collect(keys(b.faces_half_edges)), by = x -> x.id)
        if !haskey(faces_map_b, f)
            faces_map_b[f] = Face(n_face)
            n_face += 1
        end
    end

    return faces_map_a, faces_map_b
end

function update_half_edges_faces!(half_edges_faces::Dict{HalfEdge, Face}, a::PlanarGraph, Ia::Vector{Vertex}, vertices_map_a::Dict{Vertex, Vertex}, faces_map_a::Dict{Face, Face})
    for (he, f) in a.half_edges_faces
        if !(he.src in Ia && he.dst in Ia)
            new_he = HalfEdge(vertices_map_a[he.src], vertices_map_a[he.dst])
            half_edges_faces[new_he] = faces_map_a[f]
        end
    end
end

@inline function contract(a::PlanarGraph, b::PlanarGraph, Ia::Vector{Vertex}, Ib::Vector{Vertex})
    @boundscheck Ia, Ib = check_indices(a, b, Ia, Ib)
    vertices_map_a, vertices_map_b = vertices_map(a, b, Ia, Ib) # original -> new
    faces_map_a, faces_map_b = faces_map(a, b, Ia, Ib)

    half_edges_faces = Dict{HalfEdge, Face}()
    update_half_edges_faces!(half_edges_faces, a, Ia, vertices_map_a, faces_map_a)
    update_half_edges_faces!(half_edges_faces, b, Ib, vertices_map_b, faces_map_b)

    return PlanarGraph(half_edges_faces)
end

function contract(a::QuonGraph, b::QuonGraph, Ia::Vector{Vertex}, Ib::Vector{Vertex})
    Ia, Ib = check_indices(a.pg, b.pg, Ia, Ib)
    vmap_a, vmap_b = vertices_map(a.pg, b.pg, Ia, Ib)
    fmap_a, fmap_b = faces_map(a.pg, b.pg, Ia, Ib)

    new_pg = contract(a.pg, b.pg, Ia, Ib)
    new_genus = Set{Face}()
    new_pos = Dict{Vertex, Tuple{Float64, Float64}}()
    for f in a.genus
        push!(new_genus, fmap_a[f])
    end
    for f in b.genus
        push!(new_genus, fmap_b[f])
    end
    for (v, p) in a.pos
        new_pos[vmap_a[v]] = p
    end
    o_b = b.pos[Ib[2]]
    o_a = a.pos[Ia[2]]
    for (v, p) in b.pos
        if !(v in Ib)
            new_pos[vmap_b[v]] = (p[1] - o_b[1] + o_a[1], p[2] - o_b[2] + o_a[2])
        end
    end
    new_pos[vmap_a[Ia[1]]] = ((a.pos[Ia[1]][1] + b.pos[Ib[1]][1] - o_b[1] + o_a[1])/2, 
                                (a.pos[Ia[1]][2] + b.pos[Ib[1]][2] - o_b[2] + o_a[2])/2)
    new_pos[vmap_a[Ia[end]]] = ((a.pos[Ia[end]][1] + b.pos[Ib[end]][1] - o_b[1] + o_a[1])/2, 
                                (a.pos[Ia[end]][2] + b.pos[Ib[end]][2] - o_b[2] + o_a[2])/2)
    return QuonGraph(new_pg, new_pos, new_genus)
end
