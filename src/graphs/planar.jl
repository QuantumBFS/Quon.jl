"""
    planar_rz()

Generate the following planar graph for Rz:

        1 ------- 2 ----- 3 ----- 4
       / \\\\ (1)  //      ||     ||  \\
      /   \\\\    //       ||     ||   \\
     /     \\\\  //  (2)   || (3) ||    \\
    9  (6)   10          11     12 (4) 13
     \\     //  \\\\        ||     ||    /
      \\   //    \\\\       ||     ||   / 
       \\ // (5)  \\\\      ||     ||  /
        8 ------- 7 ----- 6 ----- 5
"""
function planar_rz()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 10),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 13),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 8),
        Vertex(10) => HalfEdge(10, 2),
        Vertex(11) => HalfEdge(11, 3),
        Vertex(12) => HalfEdge(12, 4),
        Vertex(13) => HalfEdge(13, 4),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(1) => HalfEdge(1, 10),
        Face(2) => HalfEdge(2, 10),
        Face(3) => HalfEdge(3, 11),
        Face(4) => HalfEdge(4, 12),
        Face(5) => HalfEdge(8, 7),
        Face(6) => HalfEdge(1, 9),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(2, 1) => Face(1),
        HalfEdge(1, 10) => Face(1),
        HalfEdge(10, 2) => Face(1),

        HalfEdge(2, 10) => Face(2),
        HalfEdge(10, 7) => Face(2),
        HalfEdge(7, 6) => Face(2),
        HalfEdge(6, 11) => Face(2),
        HalfEdge(11, 3) => Face(2),
        HalfEdge(3, 2) => Face(2),

        HalfEdge(3, 11) => Face(3),
        HalfEdge(11, 6) => Face(3),
        HalfEdge(6, 5) => Face(3),
        HalfEdge(5, 12) => Face(3),
        HalfEdge(12, 4) => Face(3),
        HalfEdge(4, 3) => Face(3),

        HalfEdge(4, 12) => Face(4),
        HalfEdge(12, 5) => Face(4),
        HalfEdge(5, 13) => Face(4),
        HalfEdge(13, 4) => Face(4),

        HalfEdge(8, 7) => Face(5),
        HalfEdge(7, 10) => Face(5),
        HalfEdge(10, 8) => Face(5),

        HalfEdge(1, 9) => Face(6),
        HalfEdge(9, 8) => Face(6),
        HalfEdge(8, 10) => Face(6),
        HalfEdge(10, 1) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end

"""
    planar_rx()

Generate the following planar graph for Rx:

         1 ----- 2 ------ 3 ----- 4
       /  ||     \\\\ (3)  //     ||  \\
      /   ||      \\\\    //      ||   \\
     /    ||  (2)  \\\\  //  (4)  ||    \\
    9 (1) 10         11         12 (5) 13
     \\    ||       //  \\\\       ||    /
      \\   ||      //    \\\\      ||   / 
       \\  ||     // (6)  \\\\     ||  /
         8 ----- 7 ------ 6 ----- 5
"""
function planar_rx()
    vertices = Dict{Vertex, HalfEdge}(
        Vertex(1) => HalfEdge(1, 10),
        Vertex(2) => HalfEdge(2, 1),
        Vertex(3) => HalfEdge(3, 2),
        Vertex(4) => HalfEdge(4, 3),
        Vertex(5) => HalfEdge(5, 13),
        Vertex(6) => HalfEdge(6, 5),
        Vertex(7) => HalfEdge(7, 6),
        Vertex(8) => HalfEdge(8, 7),
        Vertex(9) => HalfEdge(9, 8),
        Vertex(10) => HalfEdge(10, 1),
        Vertex(11) => HalfEdge(11, 2),
        Vertex(12) => HalfEdge(12, 4),
        Vertex(13) => HalfEdge(13, 4),
    )

    faces_half_edges = Dict{Face, HalfEdge}(
        Face(1) => HalfEdge(1, 9),
        Face(2) => HalfEdge(2, 1),
        Face(3) => HalfEdge(3, 2),
        Face(4) => HalfEdge(4, 3),
        Face(5) => HalfEdge(5, 13),
        Face(6) => HalfEdge(6, 11),
    )

    half_edges_faces = Dict{HalfEdge, Face}(
        HalfEdge(1, 9) => Face(1),
        HalfEdge(9, 8) => Face(1),
        HalfEdge(8, 10) => Face(1),
        HalfEdge(10, 1) => Face(1),

        HalfEdge(2, 1) => Face(2),
        HalfEdge(1, 10) => Face(2),
        HalfEdge(10, 8) => Face(2),
        HalfEdge(8, 7) => Face(2),
        HalfEdge(7, 11) => Face(2),
        HalfEdge(11, 2) => Face(2),

        HalfEdge(2, 11) => Face(3),
        HalfEdge(11, 3) => Face(3),
        HalfEdge(3, 2) => Face(3),

        HalfEdge(4, 3) => Face(4),
        HalfEdge(3, 11) => Face(4),
        HalfEdge(11, 6) => Face(4),
        HalfEdge(6, 5) => Face(4),
        HalfEdge(5, 12) => Face(4),
        HalfEdge(12, 4) => Face(4),

        HalfEdge(5, 13) => Face(5),
        HalfEdge(13, 4) => Face(5),
        HalfEdge(4, 12) => Face(5),
        HalfEdge(12, 5) => Face(5),

        HalfEdge(6, 11) => Face(6),
        HalfEdge(11, 7) => Face(6),
        HalfEdge(7, 6) => Face(6),
    )

    return PlanarGraph(vertices, faces_half_edges, half_edges_faces)
end
