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

function tait_rx(θ::T) where T
    g = planar_rx()
    p = Phase(θ, true)
    phases = Dict(7 => p, 8 => p)
    inputs = Int[1]
    outputs = Int[6]
    genuses = Set([2, 5])
    locations = Dict(
        1 => (1.0, 0.0),
        2 => (0.0, 1.0),
        3 => (1.0, 1.0),
        4 => (1.0, 2.0),
        5 => (2.0, 1.0),
        6 => (1.0, 3.0),
    )

    return Tait{Phase{T}}(g, phases, inputs, outputs, genuses, locations)
end

function tait_rz(θ::T) where T
    g = planar_rz()
    p = Phase(θ, false)
    phases = Dict(7 => p, 8 => p)
    inputs = Int[1]
    outputs = Int[5]
    genuses = Set([2, 4])
    locations = Dict(
        1 => (1.0, 0.0),
        2 => (0.0, 1.0),
        3 => (1.0, 1.0),
        4 => (2.0, 1.0),
        5 => (1.0, 2.0),
    )

    return Tait{Phase{T}}(g, phases, inputs, outputs, genuses, locations)
end

function planar_id()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 7,
        3 => 4,
        4 => 6,
        5 => 12,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 2),
        2 => HalfEdge(2, 1),
        3 => HalfEdge(1, 3),
        4 => HalfEdge(3, 1),
        5 => HalfEdge(1, 4),
        6 => HalfEdge(4, 1),
        7 => HalfEdge(2, 5),
        8 => HalfEdge(5, 2),
        9 => HalfEdge(3, 5),
        10 => HalfEdge(5, 3),
        11 => HalfEdge(4, 5),
        12 => HalfEdge(5, 4),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 3,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        11 => 0,
        8 => 0,

        1 => 1,
        7 => 1,
        10 => 1,
        4 => 1,

        3 => 2,
        9 => 2,
        12 => 2,
        6 => 2,
    )    
    next = Dict{Int, Int}(
        2 => 5,
        5 => 11,
        11 => 8,
        8 => 2,

        1 => 7,
        7 => 10,
        10 => 4,
        4 => 1,

        3 => 9,
        9 => 12,
        12 => 6,
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
    )
    id = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin,
        5, 12, 2
    )
    return id
end

function tait_id()
    g = planar_id()
    phases = Dict{Int, Phase{ComplexF64}}()
    inputs = Int[1]
    outputs = Int[5]
    genuses = Set([2, 4])
    locations = Dict(
        1 => (1.0, 0.0),
        2 => (0.0, 1.0),
        3 => (1.0, 1.0),
        4 => (2.0, 1.0),
        5 => (1.0, 2.0),
    )

    return Tait{Phase{ComplexF64}}(g, phases, inputs, outputs, genuses, locations)
end

function tait_hadamard()
    z1 = tait_rz(im*pi/2)
    x = tait_rx(im*pi/2)
    z2 = tait_rz(im*pi/2)
    contract!(z1, x)
    contract!(z1, z2)
    return z1
end

function planar_copy()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 7,
        3 => 13,
        4 => 8,
        5 => 6,
        6 => 2,
        7 => 4,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 6),
        2 => HalfEdge(6, 1),
        3 => HalfEdge(1, 7),
        4 => HalfEdge(7, 1),
        5 => HalfEdge(1, 5),
        6 => HalfEdge(5, 1),
        7 => HalfEdge(2, 4),
        8 => HalfEdge(4, 2),
        9 => HalfEdge(2, 7),
        10 => HalfEdge(7, 2),
        11 => HalfEdge(2, 6),
        12 => HalfEdge(6, 2),
        13 => HalfEdge(3, 5),
        14 => HalfEdge(5, 3),
        15 => HalfEdge(3, 7),
        16 => HalfEdge(7, 3),
        17 => HalfEdge(3, 4),
        18 => HalfEdge(4, 3),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 7,
        3 => 3,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        14 => 0,
        17 => 0,
        8 => 0,
        11 => 0,

        1 => 1,
        12 => 1,
        9 => 1,
        4 => 1,

        7 => 2,
        18 => 2,
        15 => 2,
        10 => 2,

        3 => 3,
        16 => 3,
        13 => 3,
        6 => 3,
    )    
    next = Dict{Int, Int}(
        2 => 5,
        5 => 14,
        14 => 17,
        17 => 8,
        8 => 11,
        11 => 2,

        1 => 12,
        12 => 9,
        9 => 4,
        4 => 1,

        7 => 18,
        18 => 15,
        15 => 10,
        10 => 7,

        3 => 16,
        16 => 13,
        13 => 6,
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
        15 => 16,
        16 => 15,
        17 => 18,
        18 => 17,
    )
    copy_tensor = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin,
        7, 18, 3
    )
    return copy_tensor
end

function tait_copy()
    g = planar_copy()
    phases = Dict{Int, Phase{ComplexF64}}()
    inputs = Int[1]
    outputs = Int[2, 3]
    genuses = Set([4, 5, 6])
    locations = Dict(
        1 => (1.0, 0.0),
        2 => (0.0, 2.0),
        3 => (2.0, 2.0),
        4 => (1.0, 2.0),
        5 => (2.0, 1.0),
        6 => (0.0, 1.0),
        7 => (1.0, 1.0),
    )

    return Tait{Phase{ComplexF64}}(g, phases, inputs, outputs, genuses, locations)
end