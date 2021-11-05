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
    vs_isolated = Dict{Int, Int}()
    rx = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated,
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
    vs_isolated = Dict{Int, Int}()
    rz = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated,
        5, 14, 3
    )
    return rz
end

function tait_rx(θ::T) where T
    g = planar_rx()
    p = QuonParam(θ, true)
    quon_params = Dict(7 => p, 8 => p)
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

    return Tait{QuonParam{T}}(g, quon_params, inputs, outputs, genuses, locations)
end

function tait_rz(θ::T) where T
    g = planar_rz()
    p = QuonParam(θ, false)
    quon_params = Dict(7 => p, 8 => p)
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

    return Tait{QuonParam{T}}(g, quon_params, inputs, outputs, genuses, locations)
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
    vs_isolated = Dict{Int, Int}()
    id = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated, 
        5, 12, 2
    )
    return id
end

function tait_id()
    g = planar_id()
    quon_params = Dict{Int, QuonParam{ComplexF64}}()
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

    return Tait{QuonParam{ComplexF64}}(g, quon_params, inputs, outputs, genuses, locations)
end

function tait_hadamard()
    z1 = tait_rz(im*pi/2)
    x = tait_rx(im*pi/2)
    z2 = tait_rz(im*pi/2)
    contract!(z1, x)
    contract!(z1, z2)
    z1.locations[13] = (0.0, 1.5)
    z1.locations[15] = (2.0, 1.5)
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
    vs_isolated = Dict{Int, Int}()
    copy_tensor = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated,
        7, 18, 3
    )
    return copy_tensor
end

function tait_copy()
    g = planar_copy()
    quon_params = Dict{Int, QuonParam{ComplexF64}}()
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

    return Tait{QuonParam{ComplexF64}}(g, quon_params, inputs, outputs, genuses, locations)
end

function tait_cz()
    h = tait_hadamard()
    c1 = tait_copy()
    c2 = tait_copy()
    contract!(c1, h, [3], [1])
    contract!(c1, c2, [23], [2])
    c1.locations[2] = (1.0, 3.0)
    c1.locations[6] = (0.0, 1.0)
    c1.locations[29] = (2.0, 0.0)
    c1.locations[24] = (3.0, 0.0)
    c1.locations[26] = (3.0, 3.0)
    c1.locations[28] = (4.0, 1.0)
    c1.locations[30] = (3.0, 1.0)
    c1.locations[27] = (2.0, 2.0)

    return c1
end

function planar_zero_state()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 2,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 2),
        2 => HalfEdge(2, 1),
        3 => HalfEdge(1, 2),
        4 => HalfEdge(2, 1),
        5 => HalfEdge(1, 2),
        6 => HalfEdge(2, 1),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 3,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        
        1 => 1,
        4 => 1,

        3 => 2,
        6 => 2,
    )    
    next = Dict{Int, Int}(
        2 => 5,
        5 => 2,

        1 => 4,
        4 => 1,

        3 => 6,
        6 => 3,
    )
    twin = Dict{Int, Int}(
        1 => 2,
        2 => 1,
        3 => 4,
        4 => 3,
        5 => 6,
        6 => 5,
    )
    vs_isolated = Dict{Int, Int}()
    zero_state = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated,
        2, 6, 2
    )
    return zero_state
end

function tait_zero_state()
    g = planar_zero_state()
    quon_params = Dict{Int, QuonParam{ComplexF64}}()
    inputs = Int[]
    outputs = Int[2]
    genuses = Set([1])
    locations = Dict(
        1 => (1.0, 0.0),
        2 => (1.0, 1.0),
    )

    return Tait{QuonParam{ComplexF64}}(g, quon_params, inputs, outputs, genuses, locations)
end

function planar_swap()
    v2he = Dict{Int, Int}(
        1 => 1,
        2 => 7,
        3 => 13,
        4 => 19,
        5 => 2,
        6 => 8,
        7 => 14,
        8 => 6,
        9 => 4,
        10 => 10,
        11 => 16,
        12 => 22,
        13 => 25,
        14 => 33,
        15 => 41,
        16 => 49,
        17 => 30,
    )
    half_edges = Dict{Int, HalfEdge}(
        1 => HalfEdge(1, 5),
        2 => HalfEdge(5, 1),
        3 => HalfEdge(1, 9),
        4 => HalfEdge(9, 1),
        5 => HalfEdge(1, 8),
        6 => HalfEdge(8, 1),
        7 => HalfEdge(2, 6),
        8 => HalfEdge(6, 2),
        9 => HalfEdge(2, 10),
        10 => HalfEdge(10, 2),
        11 => HalfEdge(2, 5),
        12 => HalfEdge(5, 2),
        13 => HalfEdge(3, 7),
        14 => HalfEdge(7, 3),
        15 => HalfEdge(3, 11),
        16 => HalfEdge(11, 3),
        17 => HalfEdge(3, 6),
        18 => HalfEdge(6, 3),
        19 => HalfEdge(4, 8),
        20 => HalfEdge(8, 4),
        21 => HalfEdge(4, 12),
        22 => HalfEdge(12, 4),
        23 => HalfEdge(4, 7),
        24 => HalfEdge(7, 4),
        25 => HalfEdge(13, 5),
        26 => HalfEdge(5, 13),
        27 => HalfEdge(13, 10),
        28 => HalfEdge(10, 13),
        29 => HalfEdge(13, 17),
        30 => HalfEdge(17, 13),
        31 => HalfEdge(13, 9),
        32 => HalfEdge(9, 13),
        33 => HalfEdge(14, 6),
        34 => HalfEdge(6, 14),
        35 => HalfEdge(14, 11),
        36 => HalfEdge(11, 14),
        37 => HalfEdge(14, 17),
        38 => HalfEdge(17, 14),
        39 => HalfEdge(14, 10),
        40 => HalfEdge(10, 14),
        41 => HalfEdge(15, 7),
        42 => HalfEdge(7, 15),
        43 => HalfEdge(15, 12),
        44 => HalfEdge(12, 15),
        45 => HalfEdge(15, 17),
        46 => HalfEdge(17, 15),
        47 => HalfEdge(15, 11),
        48 => HalfEdge(11, 15),
        49 => HalfEdge(16, 8),
        50 => HalfEdge(8, 16),
        51 => HalfEdge(16, 9),
        52 => HalfEdge(9, 16),
        53 => HalfEdge(16, 17),
        54 => HalfEdge(17, 16),
        55 => HalfEdge(16, 12),
        56 => HalfEdge(12, 16),
    )
    
    f2he = Dict{Int, Int}(
        0 => 2,
        1 => 1,
        2 => 3,
        3 => 7,
        4 => 9,
        5 => 13,
        6 => 15,
        7 => 19,
        8 => 21,
        9 => 29,
        10 => 27,
        11 => 35,
        12 => 43,
    )
    he2f = Dict{Int, Int}(
        2 => 0,
        5 => 0,
        20 => 0,
        23 => 0,
        14 => 0,
        17 => 0,
        8 => 0,
        11 => 0,

        1 => 1,
        26 => 1,
        31 => 1,
        4 => 1,

        3 => 2,
        52 => 2,
        49 => 2,
        6 => 2,

        7 => 3,
        34 => 3,
        39 => 3,
        10 => 3,

        9 => 4,
        28 => 4,
        25 => 4,
        12 => 4,

        13 => 5,
        42 => 5,
        47 => 5,
        16 => 5,

        15 => 6,
        36 => 6,
        33 => 6,
        18 => 6,

        19 => 7,
        50 => 7,
        55 => 7,
        22 => 7,

        21 => 8,
        44 => 8,
        41 => 8,
        24 => 8,

        29 => 9,
        54 => 9,
        51 => 9,
        32 => 9,

        27 => 10,
        40 => 10,
        37 => 10,
        30 => 10,

        35 => 11,
        48 => 11,
        45 => 11,
        38 => 11,

        43 => 12,
        56 => 12,
        53 => 12,
        46 => 12,
    )
    next = Dict{Int, Int}(
        2 => 5,
        5 => 20,
        20 => 23,
        23 => 14,
        14 => 17,
        17 => 8,
        8 => 11,
        11 => 2,

        1 => 26,
        26 => 31,
        31 => 4,
        4 => 1,

        3 => 52,
        52 => 49,
        49 => 6,
        6 => 3,

        7 => 34,
        34 => 39,
        39 => 10,
        10 => 7,

        9 => 28,
        28 => 25,
        25 => 12,
        12 => 9,

        13 => 42,
        42 => 47,
        47 => 16,
        16 => 13,

        15 => 36,
        36 => 33,
        33 => 18,
        18 => 15,

        19 => 50,
        50 => 55,
        55 => 22,
        22 => 19,

        21 => 44,
        44 => 41,
        41 => 24,
        24 => 21,

        29 => 54,
        54 => 51,
        51 => 32,
        32 => 29,

        27 => 40,
        40 => 37,
        37 => 30,
        30 => 27,

        35 => 48,
        48 => 45,
        45 => 38,
        38 => 35,

        43 => 56,
        56 => 53,
        53 => 46,
        46 => 43,
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
        19 => 20,
        20 => 19,
        21 => 22,
        22 => 21,
        23 => 24,
        24 => 23,
        25 => 26,
        26 => 25,
        27 => 28,
        28 => 27,
        29 => 30,
        30 => 29,
        31 => 32,
        32 => 31,
        33 => 34,
        34 => 33,
        35 => 36,
        36 => 35,
        37 => 38,
        38 => 37,
        39 => 40,
        40 => 39,
        41 => 42,
        42 => 41,
        43 => 44,
        44 => 43,
        45 => 46,
        46 => 45,
        47 => 48,
        48 => 47,
        49 => 50,
        50 => 49,
        51 => 52,
        52 => 51,
        53 => 54,
        54 => 53,
        55 => 56,
        56 => 55,
    )
    vs_isolated = Dict{Int, Int}()
    swap = PlanarMultigraph(v2he, half_edges,
        f2he, he2f,
        next, twin, vs_isolated, 
        17, 56, 12
    )
    return swap
end

function tait_swap()
    g = planar_swap()
    p_para = QuonParam(π/2*im, true)
    p_prop = QuonParam(π/2*im, false)
    quon_params = Dict{Int, QuonParam{ComplexF64}}(
        25 => p_para,
        26 => p_para,
        27 => p_prop,
        28 => p_prop,
        29 => p_para,
        30 => p_para,
        31 => p_prop,
        32 => p_prop,

        33 => p_prop,
        34 => p_prop,
        35 => p_para,
        36 => p_para,
        37 => p_prop,
        38 => p_prop,
        39 => p_para,
        40 => p_para,

        41 => p_para,
        42 => p_para,
        43 => p_prop,
        44 => p_prop,
        45 => p_para,
        46 => p_para,
        47 => p_prop,
        48 => p_prop,

        49 => p_prop,
        50 => p_prop,
        51 => p_para,
        52 => p_para,
        53 => p_prop,
        54 => p_prop,
        55 => p_para,
        56 => p_para,
    )
    inputs = Int[1, 4]
    outputs = Int[2, 3]
    genuses = Set([5, 6, 7, 8])
    locations = Dict(
        1 => (0.0, 0.0),
        8 => (2.0, 0.0),
        4 => (4.0, 0.0),
        9 => (1.0, 1.0),
        16 => (2.0, 1.0),
        12 => (3.0, 1.0),
        5 => (0.0, 2.0),
        13 => (1.0, 2.0),
        17 => (2.0, 2.0),
        15 => (3.0, 2.0),
        7 => (4.0, 2.0),
        10 => (1.0, 3.0),
        14 => (2.0, 3.0),
        11 => (3.0, 3.0),
        2 => (0.0, 4.0),
        6 => (2.0, 4.0),
        3 => (4.0, 4.0),
    )

    return Tait{QuonParam{ComplexF64}}(g, quon_params, inputs, outputs, genuses, locations)
end

function tait_cnot(ctrl, loc)
    z_copy = tait_copy()
    x_copy = contract!(tait_hadamard(), tait_copy())
    h = tait_hadamard()
    for v in copy(x_copy.outputs)
        contract!(x_copy, h, [v], copy(h.inputs))
    end

    if ctrl < loc
        push!(x_copy.inputs, pop!(x_copy.outputs))
        left, right = z_copy, x_copy
    elseif ctrl > loc
        push!(z_copy.inputs, pop!(z_copy.outputs))
        left, right = x_copy, z_copy
    else
        error("`ctrl` and `loc` should be different")
    end

    s = tait_swap()
    for _ in 1:(abs(ctrl - loc) - 1)
        contract!(left, s, [left.outputs[end]], [s.inputs[1]])
    end

    return contract!(left, right, [left.outputs[end]], [right.inputs[1]])
end

function tait_cz(ctrl, loc)
    @assert ctrl != loc
    left = tait_copy()
    right = tait_copy()
    h = tait_hadamard()
    contract!(left, h, [left.outputs[end]], copy(h.inputs))
    push!(right.inputs, pop!(right.outputs))
    s = tait_swap()
    for _ in 1:(abs(ctrl - loc) - 1)
        contract!(left, s, [left.outputs[end]], [s.inputs[1]])
    end

    return contract!(left, right, [left.outputs[end]], [right.inputs[1]])
end
