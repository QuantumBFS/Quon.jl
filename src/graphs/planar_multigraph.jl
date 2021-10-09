struct HalfEdge
    src::Int
    dst::Int
end
src(he::HalfEdge) = he.src
dst(he::HalfEdge) = he.dst


# NOTE: 
# isolated vertices can only appear without any other connection.
# thus we can just use a vertice -> face map for isolated vertices.
mutable struct PlanarMultigraph
    v2he::Dict{Int, Int}    # v_id -> he_id
    half_edges::Dict{Int, HalfEdge} # he_id -> he
    
    f2he::Dict{Int, Int}    # f_id -> he_id
    he2f::Dict{Int, Int}    # he_id -> f_id
    
    next::Dict{Int, Int}    # he_id -> he_id
    twin::Dict{Int, Int}    # he_id -> he_id

    vs_isolated::Dict{Int, Int} # v_id -> f_id
    
    v_max::Int
    he_max::Int
    f_max::Int
end

vertices(g::PlanarMultigraph) = vcat(collect(keys(g.v2he)), collect(keys(g.vs_isolated)))
isolated_vertices(g::PlanarMultigraph) = collect(keys(g.vs_isolated))
is_isolated(g::PlanarMultigraph, v::Integer) = haskey(g.vs_isolated, v)

faces(g::PlanarMultigraph) = sort!(collect(keys(g.f2he)))
half_edges(g::PlanarMultigraph) = sort!(collect(keys(g.half_edges)))

src(g::PlanarMultigraph, he_id::Integer) = src(g.half_edges[he_id])
dst(g::PlanarMultigraph, he_id::Integer) = dst(g.half_edges[he_id])
half_edge(g::PlanarMultigraph, he_id::Integer) = g.half_edges[he_id]
face(g::PlanarMultigraph, he_id::Integer) = g.he2f[he_id]

α(g::PlanarMultigraph, he::Integer) = g.next[he]
next(g::PlanarMultigraph, he::Integer) = α(g, he)
function prev(g::PlanarMultigraph, he::Integer)
    current_he = next(g, he)
    next_he = next(g, current_he)
    while next_he != he
        current_he = next_he
        next_he = next(g, next_he)
    end
    return current_he
end

# NOTE: rename ϕ and σ
ϕ(g::PlanarMultigraph, he::Integer) = g.twin[he]
twin(g::PlanarMultigraph, he::Integer) = ϕ(g, he)

σ(g::PlanarMultigraph, he::Integer) = twin(g, prev(g, he))
σ_inv(g::PlanarMultigraph, he::Integer) = next(g, twin(g, he))

nv(g::PlanarMultigraph) = length(g.v2he)
nf(g::PlanarMultigraph) = length(g.f2he)
nhe(g::PlanarMultigraph) = length(g.half_edges)
ne(g::PlanarMultigraph) = nhe(g) ÷ 2

function out_half_edge(g::PlanarMultigraph, v::Integer)
    is_isolated(g, v) && return 0
    return g.v2he[v]
end
surrounding_half_edge(g::PlanarMultigraph, f::Integer) = g.f2he[f]

function trace_orbit(f::Function, a::T; rev::Bool = false) where T
    next = f(a)
    perm = T[a]
    while next != a
        if rev 
            pushfirst!(perm, next)
        else
            push!(perm, next)
        end
        next = f(next)
    end
    return perm
end
function trace_face(g::PlanarMultigraph, f::Integer; safe_trace = false) 
    !safe_trace && return trace_orbit(h -> g.next[h], surrounding_half_edge(g, f))
    hes_f = Int[]
    for (he, f_he) in g.he2f
        f_he == f && push!(hes_f, he)
    end
    return hes_f
end
function trace_vertex(g::PlanarMultigraph, v::Integer)
    is_isolated(g, v) && return Int[]
    return trace_orbit(h -> σ_inv(g, h), out_half_edge(g, v); rev = true)
end
neighbors(g::PlanarMultigraph, v::Integer) = [dst(g, he) for he in trace_vertex(g, v)]

is_boundary(g::PlanarMultigraph, he_id::Integer) = (face(g, he_id) == 0)

function rem_vertex!(g::PlanarMultigraph, v::Integer; update::Bool = true)
    for he in trace_vertex(g, v)
        rem_edge!(g, he; update = update)
    end
    delete!(g.v2he, v)
    delete!(g.vs_isolated, v)
    return g
end

function rem_edge!(g::PlanarMultigraph, he_id::Integer; update::Bool = true)
    # make sure the face of he_id is an inner face
    if is_boundary(g, he_id)
        he_id = twin(g, he_id)
    end

    # handle self-loop
    if next(g, twin(g, he_id)) == twin(g, he_id)
        if next(g, he_id) == he_id  # isolated self-loop
            v_loop = src(g, he_id)
            twin_id = twin(g, he_id)
            he_in = (surrounding_half_edge(g, face(g, he_id)) == he_id) ? he_id : twin_id
            he_out = twin(g, he_in)
            f_in = face(g, he_in)
            f_out = face(g, he_out)
            for (v_iso, f_iso) in g.vs_isolated
                f_iso == f_in && (g.vs_isolated[v_iso] = f_out)
            end
            hes_f_he = trace_face(g, f_in; safe_trace = true)
            rem_face!(g, f_in)
            for he in hes_f_he
                g.he2f[he] = f_out
            end
            for he_rm in (he_in, he_out)
                delete!(g.he2f, he_rm)
                delete!(g.half_edges, he_rm)
                delete!(g.next, he_rm)
                delete!(g.twin, he_rm)
            end
            if update
                delete!(g.v2he, v_loop)
                g.vs_isolated[v_loop] = f_out
            end
            return g
        end
        he_id = twin(g, he_id)
        is_boundary(g, he_id) && (he_id = twin(g, he_id))
    end

    twin_id = twin(g, he_id)

    if update
        he_next = next(g, he_id)
        he_prev = prev(g, he_id)
        twin_next = next(g, twin_id)
        twin_prev = prev(g, twin_id)

        face_he = face(g, he_id)
        face_twin = face(g, twin_id)

        # remove a inner face
        if face_he != face_twin 
            for (v_iso, f_iso) in g.vs_isolated
                f_iso == face_he && (g.vs_isolated[v_iso] = face_twin)
            end
            hes_f_he = trace_face(g, face_he; safe_trace = true)
            rem_face!(g, face_he)
            for he in hes_f_he
                g.he2f[he] = face_twin
            end
        end

        # update f2he
        if surrounding_half_edge(g, face_twin) in (he_id, twin_id)
            new_he = nothing
            for nhe in (he_next, he_prev, twin_next, twin_prev)
                if !(nhe in (he_id, twin_id)) 
                    g.f2he[face_twin] = nhe
                    new_he = nhe
                    break
                end
            end
            new_he === nothing && error("surrounding half edge not founded")
        end

        # update v2he
        (out_half_edge(g, src(g, he_id)) == he_id) && (g.v2he[src(g, he_id)] = twin(g, he_prev))
        if out_half_edge(g, src(g, twin_id)) in (he_id, twin_id)
            if twin_id == twin_next
                g.v2he[src(g, twin_id)] = he_next
            else
                g.v2he[src(g, twin_id)] = twin(g, twin_prev)
            end
        end

        if he_next == he_id # he_id is the inner half edge of a self-loop
            g.next[twin_prev] = twin_next
        elseif twin_next == twin_id
            g.next[he_prev] = he_next
        else
            g.next[he_prev] = twin_next
            g.next[twin_prev] = he_next
        end

        if he_next == twin_id
            g.vs_isolated[src(g, he_next)] = face(g, he_next)
            delete!(g.v2he, src(g, he_next))
        end
        if twin_next == he_id
            g.vs_isolated[src(g, twin_next)] = face(g, twin_next)
            delete!(g.v2he, src(g, twin_next))
        end
    end
    delete!(g.next, he_id)
    delete!(g.next, twin_id)
    delete!(g.half_edges, he_id)
    delete!(g.half_edges, twin_id)
    delete!(g.twin, he_id)
    delete!(g.twin, twin_id)
    delete!(g.he2f, he_id)
    delete!(g.he2f, twin_id)

    return g
end

function rem_face!(g::PlanarMultigraph, f::Integer)
    f == 0 && error("Face 0 is for the boundary. It can not be removed.")
    half_edges_f = trace_face(g, f)
    for he in half_edges_f
        delete!(g.he2f, he)
    end
    delete!(g.f2he, f)
    return g
end

function merge_graph!(A::PlanarMultigraph, B::PlanarMultigraph)
    for (v, he_id) in B.v2he
        A.v2he[v + A.v_max] = he_id + A.he_max
    end
    for (he_id, he) in B.half_edges
        A.half_edges[he_id + A.he_max] = HalfEdge(src(he) + A.v_max, dst(he) + A.v_max)
    end
    for (f, he_id) in B.f2he
        if f != 0
            A.f2he[f + A.f_max] = he_id + A.he_max
        end
    end
    for (he_id, f) in B.he2f
        if f != 0
            A.he2f[he_id + A.he_max] = f + A.f_max
        else
            A.he2f[he_id + A.he_max] = 0
        end
    end
    for (curr, next) in B.next
        A.next[curr + A.he_max] = next + A.he_max
    end
    for (curr, twin) in B.twin
        A.twin[curr + A.he_max] = twin + A.he_max
    end
    for (v, f) in B.vs_isolated
        A.vs_isolated[v + A.v_max] = (f == 0) ? 0 : (f + A.f_max)
    end

    A.v_max += B.v_max
    A.he_max += B.he_max
    A.f_max += B.f_max

    return A
end

function check_faces(g::PlanarMultigraph)
    for f in faces(g)
        hes_f = trace_face(g, f)
        for he in hes_f
            face(g, he) == f || return false
        end
    end
    return true
end
function check_vertices(g::PlanarMultigraph)
    for v in vertices(g)
        hes_v = trace_vertex(g, v)
        for he in hes_v
            src(g, he) == v || return false
        end
    end
    return true
end
function check_combinatorial_maps(g::PlanarMultigraph)
    for he in half_edges(g)
        (he == α(g, ϕ(g, σ(g, he)))) || return false
    end
    return true
end

function update_face!(g::PlanarMultigraph, he_id)
    face_id = face(g, he_id)
    fs_rm = Int[]
    g.f2he[face_id] = he_id
    curr_he = next(g, he_id)
    while curr_he != he_id
        if g.he2f[curr_he] != face_id
            push!(fs_rm, g.he2f[curr_he])
            g.he2f[curr_he] = face_id
        end
        curr_he = next(g, curr_he)
    end
    for (v, f) in g.vs_isolated
        (f in fs_rm) && (g.vs_isolated[v] = face_id)
    end
    return g
end

function add_vertex!(g::PlanarMultigraph, f::Integer)
    haskey(g.f2he, f) || return 0
    g.v_max += 1
    v = g.v_max
    g.vs_isolated[v] = f
    return v
end

function add_edge_isolated_1!(g::PlanarMultigraph, v1::Integer, v2::Integer, f::Integer)
    f == g.vs_isolated[v1] || return (0, 0)
    hes = trace_vertex(g, v2) 
    he2_in = 0
    he2_out = 0
    for he in hes
        if face(g, twin(g, he)) == f
            he2_in = twin(g, he)
            he2_out = next(g, he2_in)
            break
        end
    end
    he2_in * he2_out != 0 || return (0, 0)

    g.he_max += 2
    new_he1 = g.he_max - 1
    new_he2 = g.he_max
    g.v2he[v1] = new_he2
    g.twin[new_he1] = new_he2
    g.twin[new_he2] = new_he1
    g.next[he2_in] = new_he1
    g.next[new_he1] = new_he2
    g.next[new_he2] = he2_out
    g.he2f[new_he1] = f
    g.he2f[new_he2] = f
    g.half_edges[new_he1] = HalfEdge(v2, v1)
    g.half_edges[new_he2] = HalfEdge(v1, v2)
    delete!(g.vs_isolated, v1)

    return (new_he1, new_he2)
end
function add_edge_isolated_2!(g::PlanarMultigraph, v1::Integer, v2::Integer, f::Integer)
    f == g.vs_isolated[v1] == g.vs_isolated[v2] || return (0, 0)
    g.he_max += 2
    new_he1 = g.he_max - 1
    new_he2 = g.he_max
    g.twin[new_he1] = new_he2
    g.twin[new_he2] = new_he1
    g.next[new_he1] = new_he2
    g.next[new_he2] = new_he1
    g.v2he[v1] = new_he1
    g.v2he[v2] = new_he2
    g.he2f[new_he1] = f
    g.he2f[new_he2] = f
    g.half_edges[new_he1] = HalfEdge(v1, v2)
    g.half_edges[new_he2] = HalfEdge(v2, v1)
    delete!(g.vs_isolated, v1)
    delete!(g.vs_isolated, v2)

    return (new_he1, new_he2)
end

function add_edge!(g::PlanarMultigraph, v1::Integer, v2::Integer, f::Integer)
    if is_isolated(g, v1)
        if is_isolated(g, v2)
            return add_edge_isolated_2!(g, v1, v2, f)
        else
            return add_edge_isolated_1!(g, v1, v2, f)
        end
    elseif is_isolated(g, v2)
        return add_edge_isolated_1!(g, v2, v1, f)
    end
    hes_f = trace_face(g, f)
    he1_in, he1_out, he2_in, he2_out = (0,0,0,0)
    for he in hes_f
        dst(g, he) == v1 && (he1_in = he)
        src(g, he) == v1 && (he1_out = he)
        dst(g, he) == v2 && (he2_in = he)
        src(g, he) == v2 && (he2_out = he)
    end
    he1_in * he1_out * he2_in * he2_out != 0 || return (0, 0)
    new_he1 = g.he_max + 1
    new_he2 = g.he_max + 2
    g.he_max += 2
    g.twin[new_he1] = new_he2
    g.twin[new_he2] = new_he1
    g.half_edges[new_he1] = HalfEdge(v1, v2)
    g.half_edges[new_he2] = HalfEdge(v2, v1)
    g.next[he1_in] = new_he1
    g.next[new_he1] = he2_out
    g.next[he2_in] = new_he2
    g.next[new_he2] = he1_out
    g.he2f[new_he1] = f
    g.f2he[f] = new_he1
    g.f_max += 1
    g.he2f[new_he2] = g.f_max
    g.f2he[g.f_max] = new_he2
    update_face!(g, new_he2)
    return (new_he1, new_he2)
end

function contract_edge!(g::PlanarMultigraph, he_id::Integer)
    twin_id = twin(g, he_id)
    he_prev = prev(g, he_id)
    he_next = next(g, he_id)
    twin_prev = prev(g, twin_id)
    twin_next = next(g, twin_id)
    
    v1 = src(g, he_id)
    v2 = dst(g, he_id)
    if v1 == v2
        rem_edge!(g, he_id; update = true)
        return (v1, v2)
    end
    if length(trace_vertex(g, v1)) == 1
        if length(trace_vertex(g, v2)) > 1
            (v11, v22) = contract_edge!(g, twin_id)
            return (v22, v11)
        else
            # 2 isolated vertices
            error("TODO: fix this")
        end
    end
    
    # update out half edge of v1
    out_half_edge(g, v1) == he_id && (g.v2he[v1] = twin_next)
    for he in trace_vertex(g, v2)
        v0 = dst(g, he)
        g.half_edges[he] = HalfEdge(v1, v0)
        g.half_edges[twin(g, he)] = HalfEdge(v0, v1)
    end

    if he_next == twin_id
        g.next[he_prev] = twin_next
        g.f2he[face(g, he_id)] = he_prev
    else
        g.next[he_prev] = he_next
        g.next[twin_prev] = twin_next
        g.f2he[face(g, he_id)] = he_prev
        g.f2he[face(g, twin_id)] = twin_prev
    end
    delete!(g.next, he_id)
    delete!(g.next, twin_id)
    delete!(g.half_edges, he_id)
    delete!(g.half_edges, twin_id)
    delete!(g.twin, he_id)
    delete!(g.twin, twin_id)
    delete!(g.he2f, he_id)
    delete!(g.he2f, twin_id)
    delete!(g.v2he, v2)
    # v2 is removed
    return (v1, v2)
end

has_vertex(g::PlanarMultigraph, v) = haskey(g.v2he, v) || haskey(g.vs_isolated, v)
has_half_edge(g::PlanarMultigraph, he) = haskey(g.half_edges, he)
has_face(g::PlanarMultigraph, f) = haskey(g.f2he, f)