using Test
using Quon

a = planar_rz()
b = planar_rx()
Ia = Vertex[8, 7, 6, 5]
Ib = Vertex[1, 2, 3, 4]

@testset "vertices_map" begin
    map_a, map_b = vertices_map(a, b, Ia, Ib)
    for (va, vb) in zip(Ia, Ib)
        @test map_a[va] == map_b[vb]
    end
end

@testset "faces_map" begin
    map_a, map_b = faces_map(a, b, Ia, Ib)
    @test map_a[Face(5)] == map_b[Face(2)]
    @test map_a[Face(2)] == map_b[Face(3)]
    @test map_a[Face(3)] == map_b[Face(5)]
end
