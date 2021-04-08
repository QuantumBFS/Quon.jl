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

    return q
end

function rem_edge!(q::QuonTait, he_id::Integer; update::Bool = true)
    twin_id = twin(q, he_id)
    rem_edge!(q, he_id; update = update)
    delete!(q.phases, he_id)
    delete!(q.phases, twin_id)

    return q
end

function merge_graph!(A::QuonTait, B::QuonTait)
    v_max_A = A.g.v_max
    he_max_A = A.g.he_max
    f_max_A = A.g.f_max
    merge_graph!(A.g, B.g)
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

    return A
end

Base.copy(q::QuonTait) =
    QuonTait(copy(q.g), copy(q.phases), copy(q.inputs), copy(q.outputs), copy(q.genuses))

phase(q::QuonTait, he_id::Integer) = q.phases[he_id]
change_direction!(g::QuonTait, e_id::Integer) = change_direction!(g.phases[e_id])

function planar_rz()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 9,
        3 => 4,
        4 => 6,
        5 => 14,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 2),
        2 => HalfEdge(2, 1),
        3 => HalfEdge(1, 3),
        4 => HalfEdge(3, 1),
        5 => HalfEdge(1, 4),
        6 => HalfEdge(4, 1),
        7 => HalfEdge(2, 3),
        8 => HalfEdge(3, 2),
        9 => HalfEdge(2, 5),
        10 => HalfEdge(5, 2),
        11 => HalfEdge(3, 5),
        12 => HalfEdge(5, 3),
        13 => HalfEdge(4, 5),
        14 => HalfEdge(5, 4),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 8,
        3 => 3,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        13 => 0,
        10 => 0,

        1 => 1,
        7 => 1,
        4 => 1,

        8 => 2,
        9 => 2,
        12 => 2,
        
        3 => 3,
        11 => 3,
        14 => 3,
        6 => 3,
    )    
    next = Dict{Int, Int}(
        2 => 5,
        5 => 13,
        13 => 10,
        10 => 2,

        1 => 7,
        7 => 4,
        4 => 1,

        8 => 9,
        9 => 12,
        12 => 8,

        3 => 11,
        11 => 14,
        14 => 6,
        6 => 3,
    )
    twin = Dict{Int, Int}(
        1 => 2,
        2 => 1,
        3 => 4,
        4 => 3,
        5 => 6,
        6 => 5,
        7 => 8,
        8 => 7,
        9 => 10,
        10 => 9,
        11 => 12,
        12 => 11,
        13 => 14,
        14 => 13,
    )
    rz = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin,
        5, 14, 3
    )
    return rz
end

function planar_rx()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 9,
        3 => 4,
        4 => 8,
        5 => 6,
        6 => 14,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 2),
        2 => HalfEdge(2, 1),
        3 => HalfEdge(1, 3),
        4 => HalfEdge(3, 1),
        5 => HalfEdge(1, 5),
        6 => HalfEdge(5, 1),
        7 => HalfEdge(3, 4),
        8 => HalfEdge(4, 3),
        9 => HalfEdge(2, 6),
        10 => HalfEdge(6, 2),
        11 => HalfEdge(4, 6),
        12 => HalfEdge(6, 4),
        13 => HalfEdge(5, 6),
        14 => HalfEdge(6, 5),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 3,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        13 => 0,
        10 => 0,

        1 => 1,
        9 => 1,
        12 => 1,
        8 => 1,
        4 => 1,

        3 => 2,
        7 => 2,
        11 => 2,
        14 => 2,
        6 => 2,
    )    
    next = Dict{Int, Int}(
        2 => 5,
        5 => 13,
        13 => 10,
        10 => 2,

        1 => 9,
        9 => 12,
        12 => 8,
        8 => 4,
        4 => 1,

        3 => 7,
        7 => 11,
        11 => 14,
        14 => 6,
        6 => 3,
    )
    twin = Dict{Int, Int}(
        1 => 2,
        2 => 1,
        3 => 4,
        4 => 3,
        5 => 6,
        6 => 5,
        7 => 8,
        8 => 7,
        9 => 10,
        10 => 9,
        11 => 12,
        12 => 11,
        13 => 14,
        14 => 13,
    )
    rx = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin,
        6, 14, 2
    )
    return rx
end

function tait_rx(θ::T) where T
    g = planar_rx()
    p = Phase(θ, true)
    phases = Dict(7 => p, 8 => p)
    inputs = Int[1]
    outputs = Int[6]
    genuses = Set([2, 5])

    return QuonTait{Phase{T}}(g, phases, inputs, outputs, genuses)
end

function tait_rz(θ::T) where T
    g = planar_rz()
    p = Phase(θ, false)
    phases = Dict(7 => p, 8 => p)
    inputs = Int[1]
    outputs = Int[5]
    genuses = Set([2, 4])

    return QuonTait{Phase{T}}(g, phases, inputs, outputs, genuses)
end

function contract!(A::QuonTait, B::QuonTait, va::Vector{Int}, vb::Vector{Int})
    merge_graph!(A, B)
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
            !(f in fs_a) && delete!(q.g.f2he, f)
        end
        rem_vertex!(q, a; update = false)
        rem_vertex!(q, b; update = false)
    end
    return q
end

# function shift_adjlist!(adj, e_id)
#     i_a = findfirst(isequal(e_id), adj)
#     perm = [i_a+1:length(adj)..., 1:i_a...]
#     permute!(adj, perm)
#     return adj
# end
