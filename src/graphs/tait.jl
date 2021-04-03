mutable struct Phase{T <: Number}
    param::T
    isparallel::Bool
end
function change_direction!(p::Phase)
    p.param = change_direction(p.param)
    p.isparallel = !p.isparallel
    return p
end

# src < dst
struct Edge
    src::Int
    dst::Int
    function Edge(a::Integer, b::Integer)
        a <= b && return new(a, b)
        return new(b, a)
    end
end
src(e::Edge) = e.src
dst(e::Edge) = e.dst
another_end(e::Edge, v::Integer) = (v == e.src) ? e.dst : e.src

mutable struct TaitGraph{P <: Phase}
    adjlist::Dict{Int, Vector{Int}} # vertex_id => adjacent edge_ids
    edges::Dict{Int, Edge}  # edge_id => edge
    phases::Dict{Int, P} # edge_id => phase
    boundaries::Vector{Int} # list of vertex_id
    vertex_max::Int
    edge_max::Int
end

src(g::TaitGraph, e_id::Integer) = src(g.edges[e_id])
dst(g::TaitGraph, e_id::Integer) = dst(g.edges[e_id])
edge(g::TaitGraph, e_id::Integer) = g.edges[e_id]
adjacent_edges(g::TaitGraph, v_id::Integer) = g.adjlist[v_id]
neighbors(g::TaitGraph, v_id::Integer) = [another_end(edge(g, e_id), v_id) for e_id in adjacent_edges(g, v_id)]

nv(g::TaitGraph) = length(g.adjlist)
ne(g::TaitGraph) = length(g.edges)

function rem_vertex!(G::TaitGraph, v::Integer; keep_edge::Bool = false)
    if !keep_edge
        adj_v = adjacent_edges(G, v)
        nb_v = neighbors(G, v)

        for e_id in adj_v
            u = another_end(edge(G, e_id), v)
            deleteat!(G.adjlist[u], findall(isequal(e_id), G.adjlist[u]))
            delete!(G.edges, e_id)
            delete!(G.phases, e_id)
        end
    end
    delete!(G.adjlist, v)
    deleteat!(G.boundaries, findall(isequal(v), G.boundaries))

    return G
end

function rem_edge!(G::TaitGraph, e_id::Integer)
    s = src(G, e_id)
    d = dst(G, e_id)
    delete!(G.edges, e_id)
    delete!(G.phases, e_id)
    deleteat!(G.adjlist[s], findall(isequal(e_id), G.adjlist[s]))    
    deleteat!(G.adjlist[d], findall(isequal(e_id), G.adjlist[d]))    

    return G
end

Base.copy(G::TaitGraph) =
    TaitGraph(copy(G.adjlist), copy(G.edges), copy(G.phases), copy(G.boundaries), G.vertex_max, G.edge_max)

change_direction!(g::TaitGraph, e_id::Integer) = change_direction!(g.phases[e_id])

function tait_rx(θ::T) where T
    adjlist = Dict{Int, Vector{Int}}(
        1=>[1, 2],
        2=>[3, 4],
        3=>[5, 3],
        4=>[6, 7],
        5=>[2, 4, 6],
        6=>[7, 5, 1],
    )
    edges = Dict{Int, Edge}(
        1 => Edge(1, 6),
        2 => Edge(1, 5),
        3 => Edge(2, 3),
        4 => Edge(2, 5),
        5 => Edge(3, 6),
        6 => Edge(4, 5),
        7 => Edge(4, 6),
    )

    phases = Dict{Int, Phase{T}}(
        3 => Phase(θ, true),
    )
    boundaries = [5, 6]
    vmax = 6
    emax = 7
    return TaitGraph(adjlist, edges, phases, boundaries, vmax, emax)
end

function tait_rz(θ::T) where T
    adjlist = Dict{Int, Vector{Int}}(
        1=>[1, 2],
        2=>[3, 4, 5],
        3=>[7, 4, 6],
        4=>[2, 5, 7],
        5=>[6, 3, 1],
    )
    edges = Dict{Int, Edge}(
        1 => Edge(1, 5),
        2 => Edge(1, 4),
        3 => Edge(2, 5),
        4 => Edge(2, 3),
        5 => Edge(2, 4),
        6 => Edge(3, 5),
        7 => Edge(3, 4),
    )

    phases = Dict{Int, Phase{T}}(
        4 => Phase(θ, false)
    )
    boundaries = [4, 5]
    vmax = 5
    emax = 7
    return TaitGraph(adjlist, edges, phases, boundaries, vmax, emax)
end

function contract!(A::TaitGraph, B::TaitGraph, va::Vector{Int}, vb::Vector{Int})
    vb .+= A.vertex_max
    merge_graph!(A, B)
    contract_boundary_vertices!(A, va, vb)
    return A
end

merge_graph(A::TaitGraph, B::TaitGraph) = merge_graph!(copy(A), B)

function merge_graph!(A::TaitGraph, B::TaitGraph)
    for (v_id, e_ids) in B.adjlist
        A.adjlist[v_id + A.vertex_max] = e_ids .+ A.edge_max
    end

    for (e_id, e) in B.edges
        A.edges[e_id + A.edge_max] = Edge(e.src + A.vertex_max, e.dst + A.vertex_max)
    end
    for (e_id, p) in B.phases
        A.phases[e_id + A.edge_max] = p
    end

    append!(A.boundaries, B.boundaries .+ A.vertex_max)
    A.vertex_max += B.vertex_max
    A.edge_max += B.edge_max
    return A
end

# function shift_vertex!(g::TaitGraph, n::Integer)
#     n > 0 || return g
#     g.vertex_max += n
#     for (e_id, e) in edges
#         edges[e_id] = Edge(e.src + n, e.dst + n)
#     end
#     new_adjlist = Dict{Int, Vector{Int}}()
#     for v_id in sort!(collect(keys(g.adjlist)), rev = true)
#         g.adjlist[v_id + n] = g.adjlist[v_id]
#         delete!(g.adjlist, v_id)
#     end
#     g.boundaries .+= n

#     return g
# end

function contract_boundary_vertices!(G::TaitGraph, va::Vector{Int}, vb::Vector{Int})
    for (a, b) in zip(va, vb)
        adj_a = adjacent_edges(G, a)
        adj_b = adjacent_edges(G, b)
        for k in 1:length(adj_a)
            e_u = adj_a[k]
            e_v = adj_b[end-k+1]
            u = another_end(edge(G, e_u), a)
            v = another_end(edge(G, e_v), b)
            adj_u = G.adjlist[u]
            adj_v = G.adjlist[v]
            shift_adjlist!(adj_u, e_u)
            shift_adjlist!(adj_v, e_v)
            append!(adj_u, adj_v)
            
            for neighbor_v in G.adjlist[v]
                e = Edge(u, another_end(edge(G, neighbor_v), v))
                G.edges[neighbor_v] = e
            end

            rem_vertex!(G, v; keep_edge = true)
        end
        @show (a, b)
        rem_vertex!(G, a)
        rem_vertex!(G, b)
    end
    return G
end

function shift_adjlist!(adj, e_id)
    i_a = findfirst(isequal(e_id), adj)
    perm = [i_a+1:length(adj)..., 1:i_a...]
    permute!(adj, perm)
    return adj
end
