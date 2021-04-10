mutable struct Phase{T <: Number}
    param::T
    isparallel::Bool
end
Phase(p::T, isparallel::Bool) where {T <: Number} = Phase{T}(p, isparallel)

function change_direction!(p::Phase)
    p.param = change_direction(p.param)
    p.isparallel = !p.isparallel
    return p
end

mutable struct QuonTait{P <: Phase}
    g::PlanarMultigraph
    phases::Dict{Int, P}    # he_id => phase
    inputs::Vector{Int}     # list of input vertices
    outputs::Vector{Int}    # list of output vertices
    genuses::Set{Int}   # set of vertices on genuses
    locations::Dict{Int, Tuple{Float64, Float64}}   # v -> location
end

apis_1 = [:nv, :ne, :nf, :nhe, :vertices, :faces, :half_edges,
    :check_faces, :check_vertices, :check_combinatorial_maps]
apis_2 = [:src, :dst, :half_edge, :face,
    :next, :prev, :twin, :α, :ϕ, :σ, :σ_inv, :is_boundary,
    :out_half_edge, :surrounding_half_edge,
    :trace_face, :trace_vertex, :neighbors, :rem_face!
]

for api in apis_1
    @eval $api(q::QuonTait) = $api(q.g)
end
for api in apis_2
    @eval $api(q::QuonTait, id) = $api(q.g, id)
end

function rem_vertex!(q::QuonTait, v::Integer; update::Bool = true)
    rem_vertex!(q.g, v; update = update)
    deleteat!(q.inputs, findall(isequal(v), q.inputs))
    deleteat!(q.outputs, findall(isequal(v), q.outputs))
    delete!(q.genuses, v)
    delete!(q.locations, v)

    return q
end

function rem_edge!(q::QuonTait, he_id::Integer; update::Bool = true)
    twin_id = twin(q, he_id)
    rem_edge!(q, he_id; update = update)
    delete!(q.phases, he_id)
    delete!(q.phases, twin_id)

    return q
end

function merge_graph!(A::QuonTait, B::QuonTait; vertical::Bool = true, delta::Float64 = 0.0)
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

Base.copy(q::QuonTait) = QuonTait(
        copy(q.g), copy(q.phases), copy(q.inputs), copy(q.outputs),
        copy(q.genuses), copy(g.locations)
    )

phases(q::QuonTait) = q.phases
phase(q::QuonTait, he_id::Integer) = q.phases[he_id]
genuses(q::QuonTait) = sort!(collect(q.genuses))
is_genus(q::QuonTait, v::Integer) = (v in q.genuses)
change_direction!(g::QuonTait, e_id::Integer) = change_direction!(g.phases[e_id])
is_open(q::QuonTait, v::Integer) = (v in q.inputs) || (v in q.outputs)

function contract!(A::QuonTait, B::QuonTait, va::Vector{Int}, vb::Vector{Int})
    v_max_A = A.g.v_max
    merge_graph!(A, B; delta = 3.0)
    vb = vb .+ v_max_A
    contract_boundary_vertices!(A, va, vb)
    return A
end
contract!(A::QuonTait, B::QuonTait) = contract!(A, B, A.outputs, B.inputs)
merge_graph(A::QuonTait, B::QuonTait) = merge_graph!(copy(A), B)

function contract_boundary_vertices!(q::QuonTait, va::Vector{Int}, vb::Vector{Int})
    for (a, b) in zip(va, vb)
        out_a = trace_vertex(q, a)
        out_b = trace_vertex(q, b)
        ia = findfirst(he -> is_boundary(q, he), out_a)
        out_a = [out_a[ia+1:end]; out_a[1:ia]]
        ib = findfirst(he -> is_boundary(q, he), out_b)
        out_b = [out_b[ib+1:end]; out_b[1:ib]]
        reverse!(out_b)
        dst_a = [dst(q, he) for he in out_a]
        dst_b = [dst(q, he) for he in out_b]
        fs_a = [face(q, he) for he in out_a]
        fs_b = [face(q, he) for he in out_b]

        for k in 1:length(out_a)
            he_a = out_a[k]
            he_b = out_b[k]
            for he_id in trace_vertex(q, dst_b[k])
                q.g.half_edges[he_id] = HalfEdge(dst_a[k], dst(q, he_id))
                q.g.half_edges[twin(q, he_id)] = HalfEdge(dst(q, he_id), dst_a[k])
            end
            q.g.next[prev(q, twin(q, he_b))] = next(q, he_a)
            q.g.next[prev(q, twin(q, he_a))] = next(q, he_b)
        end
        for k in 1:(length(out_a) - 1)
            he_a = out_a[k]
            next_a = next(q, he_a)
            while src(q, next_a) == a || dst(q, next_a) == a
                next_a = next_a(q, he_a)
            end
            f_a = face(q, next_a)
            q.g.f2he[f_a] = next_a
            curr_he = next(q, next_a)
            while curr_he != next_a
                q.g.he2f[curr_he] = f_a
                curr_he = next(q, curr_he)
            end
        end
        for v in dst_b
            if !(v in dst_a)
                delete!(q.g.v2he, v)
                delete!(q.genuses, v)
            end
        end
        for f in fs_b
            f != 0 && delete!(q.g.f2he, f)
        end
        rem_vertex!(q, a; update = false)
        rem_vertex!(q, b; update = false)
    end
    return q
end

function right_boundary(q::QuonTait)
    v_in = q.inputs[end]
    he_in = out_half_edge(q, v_in)
    while !is_boundary(q, he_in)
        he_in = σ_inv(q, he_in)
    end
    he_bd = he_in
    right_bd = dst(q, he_bd)
    is_genus(q, right_bd) && return right_bd
end

function left_boundary(q::QuonTait)
    v_out = q.outputs[1]
    he_out = out_half_edge(q, v_out)
    while !is_boundary(q, he_out)
        he_out = σ_inv(q, he_out)
    end
    he_bd = he_out
    left_bd = dst(q, he_bd)
    is_genus(q, left_bd) && return left_bd
end

function tensor_product!(A::QuonTait, B::QuonTait; 
        right_bd_A::Integer = right_boundary(A), 
        left_bd_B::Integer = left_boundary(B))
    v_max_A = A.g.v_max
    merge_graph!(A, B; vertical = false, delta = 1.0)
    merge_boundary_vertices!(A, right_bd_A, left_bd_B + v_max_A)
    return A
end

function merge_boundary_vertices!(q::QuonTait, v1::Integer, v2::Integer)
    v1_out = out_half_edge(q, v1)
    while !is_boundary(q, v1_out)
        v1_out = σ_inv(q, v1_out)
    end
    v1_in = twin(q, σ(q, v1_out))
    v2_out = out_half_edge(q, v2)
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