mutable struct Tait{P <: Phase}
    g::PlanarMultigraph
    phases::Dict{Int, P}    # he_id => phase
    inputs::Vector{Int}     # list of input vertices
    outputs::Vector{Int}    # list of output vertices
    genuses::Set{Int}   # set of vertices on genuses
    locations::Dict{Int, Tuple{Float64, Float64}}   # v -> location
end

nv(q::Tait) = nv(q.g)
ne(q::Tait) = ne(q.g)
nf(q::Tait) = nf(q.g)
nhe(q::Tait) = nhe(q.g)
vertices(q::Tait) = vertices(q.g)
faces(q::Tait) = faces(q.g)
half_edges(q::Tait) = half_edges(q.g)
check_faces(q::Tait) = check_faces(q.g)
check_vertices(q::Tait) = check_vertices(q.g)
check_combinatorial_maps(q::Tait) = check_combinatorial_maps(q.g)
isolated_vertices(q::Tait) = isolated_vertices(q.g)

src(q::Tait, id) = src(q.g, id)
dst(q::Tait, id) = dst(q.g, id)
half_edge(q::Tait, id) = half_edge(q.g, id)
face(q::Tait, id) = face(q.g, id)
next(q::Tait, id) = next(q.g, id)
prev(q::Tait, id) = prev(q.g, id)
twin(q::Tait, id) = twin(q.g, id)
α(q::Tait, id) = α(q.g, id)
ϕ(q::Tait, id) = ϕ(q.g, id)
σ(q::Tait, id) = σ(q.g, id)
σ_inv(q::Tait, id) = σ_inv(q.g, id)
is_boundary(q::Tait, id) = is_boundary(q.g, id)
out_half_edge(q::Tait, id) = out_half_edge(q.g, id)
surrounding_half_edge(q::Tait, id) = surrounding_half_edge(q.g, id)
trace_face(q::Tait, id) = trace_face(q.g, id)
trace_vertex(q::Tait, id) = trace_vertex(q.g, id)
neighbors(q::Tait, id) = neighbors(q.g, id)
rem_face!(q::Tait, id) = rem_face!(q.g, id)
update_face!(q::Tait, id) = update_face!(q.g, id)
is_isolated(q::Tait, id) = is_isolated(q.g, id)
add_vertex!(q::Tait, id) = add_vertex!(q.g, id)

function add_edge!(q::Tait{P}, v1::Integer, v2::Integer, f::Integer, p::P) where P
    (new_he1, new_he2) = add_edge!(q.g, v1, v2, f)
    if !(0 in (new_he1, new_he2))
        q.phases[new_he1] = p
        q.phases[new_he2] = p
    end
    return (new_he1, new_he2)
end

function rem_vertex!(q::Tait, v::Integer; update::Bool = true)
    for he in trace_vertex(q, v)
        delete!(q.phases, twin(q, he))
        delete!(q.phases, he)
    end
    rem_vertex!(q.g, v; update = update)
    deleteat!(q.inputs, findall(isequal(v), q.inputs))
    deleteat!(q.outputs, findall(isequal(v), q.outputs))
    delete!(q.genuses, v)
    delete!(q.locations, v)
    
    return q
end

function rem_edge!(q::Tait, he_id::Integer; update::Bool = true)
    twin_id = twin(q, he_id)
    rem_edge!(q.g, he_id; update = update)
    delete!(q.phases, he_id)
    delete!(q.phases, twin_id)

    return q
end

function contract_edge!(q::Tait, he_id::Integer)
    twin_id = twin(q, he_id)
    (v1, v2) = contract_edge!(q.g, he_id)
    if v2 in q.genuses 
        delete!(q.genuses, v2)
        push!(q.genuses, v1)
    end
    delete!(q.phases, he_id)
    delete!(q.phases, twin_id)
    (v1 != v2) && delete!(q.locations, v2)
    return q
end

function merge_graph!(A::Tait, B::Tait; vertical::Bool = true, delta::Float64 = 0.0)
    v_max_A = A.g.v_max
    he_max_A = A.g.he_max
    f_max_A = A.g.f_max
    merge_graph!(A.g, B.g)
    ys_A = [p[2] for p in values(A.locations)]
    y_max_A = maximum(ys_A)
    xs_A = [p[1] for p in values(A.locations)]
    x_max_A = maximum(xs_A)
    for (he_id, p) in B.phases
        A.phases[he_id + he_max_A] = p
    end
    for v in B.inputs
        push!(A.inputs, v + v_max_A)
    end
    for v in B.outputs
        push!(A.outputs, v + v_max_A)
    end
    for v in B.genuses
        push!(A.genuses, v + v_max_A)
    end
    if vertical
        for (v, p) in B.locations
            A.locations[v + v_max_A] = (p[1], p[2] + y_max_A + 1 - delta)
        end
    else
        for (v, p) in B.locations
            A.locations[v + v_max_A] = (p[1] + x_max_A + 1 - delta, p[2])
        end
    end

    return A
end

Base.copy(q::Tait) = Tait(
        copy(q.g), copy(q.phases), copy(q.inputs), copy(q.outputs),
        copy(q.genuses), copy(g.locations)
    )

phases(q::Tait) = q.phases
phase(q::Tait, he_id::Integer) = q.phases[he_id]
genuses(q::Tait) = sort!(collect(q.genuses))
is_genus(q::Tait, v::Integer) = (v in q.genuses)
function change_direction!(q::Tait, e_id::Integer) 
    p = change_direction(q.phases[e_id])
    q.phases[e_id] = p
    q.phases[twin(q, e_id)] = p
    return q
end

"""
    is_open_vertex(q::Tait, v)

Returns `true` if `v` represent a special vertex for open edges.
"""
is_open_vertex(q::Tait, v::Integer) = (v in q.inputs) || (v in q.outputs)

function contract!(A::Tait, B::Tait, va::Vector{Int}, vb::Vector{Int})
    v_max_A = A.g.v_max
    merge_graph!(A, B; delta = 3.0)
    vb = vb .+ v_max_A
    contract_boundary_vertices!(A, va, vb)
    return A
end
contract!(A::Tait, B::Tait) = contract!(A, B, copy(A.outputs), copy(B.inputs))
merge_graph(A::Tait, B::Tait) = merge_graph!(copy(A), B)

function contract_boundary_vertices!(q::Tait, va::Vector{Int}, vb::Vector{Int})
    for (a, b) in zip(va, vb)
        out_a = trace_vertex(q, a)
        out_b = trace_vertex(q, b)
        ia = findfirst(he -> is_boundary(q, he), out_a)
        out_a = [out_a[ia+1:end]; out_a[1:ia]]
        ib = findfirst(he -> is_boundary(q, he), out_b)
        out_b = [out_b[ib+1:end]; out_b[1:ib]]
        reverse!(out_b)

        tp_a = [prev(q, twin(q, he)) for he in out_a]
        tp_b = [prev(q, twin(q, he)) for he in out_b]
        dst_a = [dst(q, he) for he in out_a]
        dst_b = [dst(q, he) for he in out_b]
        fs_a = [face(q, he) for he in out_a]
        fs_b = [face(q, he) for he in out_b]

        # update surrounding half edge
        ts_a = [twin(q, he) for he in out_a]
        ts_b = [twin(q, he) for he in out_b]
        ts = vcat(ts_a, ts_b)
        for he in ts
            f_he = face(q, he)
            he_f = surrounding_half_edge(q, f_he)
            while he_f in ts
                he_f = next(q, he_f)
            end
            q.g.f2he[f_he] = he_f
        end
        if dst_a[end] == dst_b[end] # this will create a self-loop
            # face(out_a[end]) should be 0
            he0 = surrounding_half_edge(q, 0)
            while (he0 in ts) || he0 == out_a[end]
                he0 = next(q, he0)
            end
            q.g.f2he[0] = he0
            q.g.f_max += 1
            new_f = q.g.f_max
            q.g.f2he[new_f] = out_a[end]
            q.g.he2f[out_a[end]] = new_f
        end

        # update faces
        for k in 1:(length(out_a)-1)
            if fs_a[k] != fs_b[k+1]
                for he in trace_face(q, fs_b[k+1])
                    q.g.he2f[he] = fs_a[k]
                end
                delete!(q.g.f2he, fs_b[k+1])
            end
        end

        # connecting 2 planar graphs
        for k in 1:length(out_a)
            t_a = ts_a[k]
            t_b = ts_b[k]
            q.g.twin[out_a[k]] = out_b[k]
            q.g.twin[out_b[k]] = out_a[k]
            q.g.next[tp_a[k]] = out_b[k]
            q.g.next[tp_b[k]] = out_a[k]
            q.g.half_edges[out_a[k]] = HalfEdge(dst_b[k], dst_a[k])
            q.g.half_edges[out_b[k]] = HalfEdge(dst_a[k], dst_b[k])
            q.g.v2he[dst_a[k]] = out_b[k]
            q.g.v2he[dst_b[k]] = out_a[k]
            delete!(q.g.he2f, t_a)
            delete!(q.g.he2f, t_b)
            delete!(q.g.next, t_a)
            delete!(q.g.next, t_b)
            delete!(q.g.twin, t_a)
            delete!(q.g.twin, t_b)
            delete!(q.g.half_edges, t_a)
            delete!(q.g.half_edges, t_b)
        end

        for v_rm in (a, b)
            delete!(q.g.v2he, v_rm)
            deleteat!(q.inputs, findall(isequal(v_rm), q.inputs))
            deleteat!(q.outputs, findall(isequal(v_rm), q.outputs))
            delete!(q.genuses, v_rm)
            delete!(q.locations, v_rm)
        end
        
        for k in 1:length(out_a)
            contract_edge!(q, out_a[k])
        end
    end
    return q
end

function right_boundary(q::Tait)
    v_in = q.inputs[end]
    he_in = out_half_edge(q, v_in)
    @assert he_in != 0
    while !is_boundary(q, he_in)
        he_in = σ_inv(q, he_in)
    end
    he_bd = he_in
    right_bd = dst(q, he_bd)
    is_genus(q, right_bd) && return right_bd
end

function left_boundary(q::Tait)
    v_out = q.outputs[1]
    he_out = out_half_edge(q, v_out)
    @assert he_out != 0
    while !is_boundary(q, he_out)
        he_out = σ_inv(q, he_out)
    end
    he_bd = he_out
    left_bd = dst(q, he_bd)
    is_genus(q, left_bd) && return left_bd
end

function tensor_product!(A::Tait, B::Tait; 
        right_bd_A::Integer = right_boundary(A), 
        left_bd_B::Integer = left_boundary(B))
    v_max_A = A.g.v_max
    merge_graph!(A, B; vertical = false, delta = 1.0)
    merge_boundary_vertices!(A, right_bd_A, left_bd_B + v_max_A)
    return A
end

function merge_boundary_vertices!(q::Tait, v1::Integer, v2::Integer)
    v1_out = out_half_edge(q, v1)
    @assert v1_out != 0
    while !is_boundary(q, v1_out)
        v1_out = σ_inv(q, v1_out)
    end
    v1_in = twin(q, σ(q, v1_out))
    v2_out = out_half_edge(q, v2)
    @assert v2_out != 0
    while !is_boundary(q, v2_out)
        v2_out = σ_inv(q, v2_out)
    end
    v2_in = twin(q, σ(q, v2_out))
    for he_id in trace_vertex(q, v2)
        twin_id = twin(q, he_id)
        he = half_edge(q, he_id)
        q.g.half_edges[he_id] = HalfEdge(v1, he.dst)
        q.g.half_edges[twin_id] = HalfEdge(he.dst, v1)
    end
    q.g.next[v1_in] = v2_out
    q.g.next[v2_in] = v1_out
    delete!(q.g.v2he, v2)
    delete!(q.genuses, v2)
    delete!(q.locations, v2)
    return q
end