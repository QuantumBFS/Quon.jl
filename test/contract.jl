using Test
using Quon

a = planar_rz()
b = planar_rx()
Ia = Vertex[9, 8, 7, 6, 5, 13]
Ib = Vertex[9, 1, 2, 3, 4, 13]

@testset "vertices_map" begin
    map_a, map_b = vertices_map(a, b, Ia, Ib)
    for (va, vb) in zip(Ia, Ib)
        @test map_a[va] == map_b[vb]
    end
end

@testset "faces_map" begin
    map_a, map_b = faces_map(a, b, Ia, Ib)
    @test map_a[Face(6)] == map_b[Face(1)]
    @test map_a[Face(5)] == map_b[Face(2)]
    @test map_a[Face(2)] == map_b[Face(3)]
    @test map_a[Face(3)] == map_b[Face(4)]
    @test map_a[Face(4)] == map_b[Face(5)]
end

contract(a, b, Ia, Ib)