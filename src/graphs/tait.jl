struct Phase{T <: Number}
    param::T
    isparallel::Bool
end

# src < dst
struct Edge
    src::Int
    dst::Int
end

mutable struct TaitGraph{P <: Phase}
    adjlist::Dict{Int, Vector{Int}} # vertex_id => adjacent edge_ids
    edges::Dict{Int, Edge}  # edge_id => edge
    phases::Dict{Int, P} # edge_id => phase
    boundaries::Vector{Int} # list of vertex_id
    vertex_max::Int
    edge_max::Int
end

Base.copy(G::TaitGraph) =
    TaitGraph(copy(G.adjlist), copy(G.edges), copy(G.phases), copy(G.boundaries), G.vertex_max, G.edge_max)

function tait_rx(θ::T) where T
    adjlist = Dict{Int, Vector{Int}}(
        1=>[6, 5],
        2=>[3, 5],
        3=>[6, 2],
        4=>[5, 6],
        5=>[1, 2, 4],
        6=>[4, 3, 1],
    )

    phases = Dict{MultiEdge, Vector{Phase{T}}}(
        MultiEdge(2, 3) => [Phase(θ, true)],
    )
    boundaries = Set([5, 6])
    idmax = 6
    return TaitGraph(adjlist, phases, boundaries, idmax)
end

function tait_rz(θ::T) where T
    adjlist = Dict{Int, Vector{Int}}(
        1=>[5, 4],
        2=>[5, 3, 4],
        3=>[4, 2, 5],
        4=>[1, 2, 3],
        5=>[3, 2, 1],
    )

    phases = Dict{MultiEdge, Vector{Phase{T}}}(
        MultiEdge(2, 3) => [Phase(θ, false)]
    )
    boundaries = Set([4, 5])
    idmax = 5
    return TaitGraph(adjlist, phases, boundaries, idmax)
end

function contract!(A::TaitGraph, B::TaitGraph, va::Vector{Int}, vb::Vector{Int})
    vb .+= A.idmax
    merge_graph!(A, B)
    contract_boundary_vertices!(A, va, vb)
    return A
end

merge_graph(A::TaitGraph, B::TaitGraph) = merge_graph!(copy(A), B)

function merge_graph!(A::TaitGraph, B::TaitGraph)
    for (vertex, adj) in B.adjlist
        A.adjlist[vertex + A.idmax] = adj
    end

    for (edge, phase) in B.phases
        e = MultiEdge(edge.src + A.idmax, edge.dst + A.idmax)
        A.phases[e] = phase
    end

    append!(A.boundaries, B.boundaries .+ A.idmax)
    A.idmax += B.idmax
    return A
end

function contract_boundary_vertices(G::TaitGraph, va::Vector{Int}, vb::Vector{Int})
    for (a, b) in zip(va, vb)
        adj_a = G.adjlist[a]
        adj_b = G.adjlist[b]
        for k in 1:length(adj_a)
            u = adj_a[k]
            v = adj_b[end-k+1]
            adj_u = G.adjlist[u]
            adj_v = G.adjlist[v]
            shift_adjlist!(adj_u, a)
            shift_adjlist!(adj_v, b)
            append!(adj_u, adj_v)
            
            for neighbor_v in G.adjlist[v]
                e = v < neighbor_v ? MultiEdge(v, neighbor_v) : MultiEdge(neighbor_v, v)
                G.phases[e]
            end

            rem_vertex!(G, v)
            delete!(G.adjlist, v)
        end
    end
end

function shift_adjlist!(adj, a)
    i_a = findfirst(isequal(a), adj)
    perm = [length(adj_u):i_a+1..., 1:i_a...]
    permute!(adj, perm)
    return adj
end
