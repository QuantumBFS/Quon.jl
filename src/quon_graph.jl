struct QuonGraph
    pg::PlanarGraph
    pos::Dict{Vertex, Tuple{Float64, Float64}}
    genus::Set{Face}
    # meta    # recording params...
    # medial  # generating medial graphs dynamically
end

function quon_rz()
    pg = planar_rz()
    pos = Dict(
        zip(Vertex[i for i = 1:13], 
            Tuple{Float64, Float64}[(1,0), (2,0), (3,0), (4,0), 
                (4,2), (3,2), (2,2), (1,2), 
                (0,1), (1.5,1), (3,1), (4,1), (5,1)]
        )
    )
    genus = Set([Face(4), Face(6)])
    return QuonGraph(pg, pos, genus)
end

function quon_rx()
    pg = planar_rx()
    pos = Dict(
        zip(Vertex[i for i = 1:13], 
            Tuple{Float64, Float64}[(1,0), (2,0), (3,0), (4,0), 
                (4,2), (3,2), (2,2), (1,2), 
                (0,1), (1,1), (2.5,1), (4,1), (5,1)]
        )
    )
    genus = Set([Face(4), Face(6)])
    return QuonGraph(pg, pos, genus)
end
