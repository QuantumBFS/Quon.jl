using Test
using Quon

a = planar_rz()
b = planar_rx()
Ia = Vertex[9, 8, 7, 6, 5, 13]
Ib = Vertex[9, 1, 2, 3, 4, 13]

@testset "vertices_map" begin
    vmap_a, vmap_b = vertices_map(a, b, Ia, Ib)
    for (va, vb) in zip(Ia, Ib)
        @test vmap_a[va] == vmap_b[vb]
    end
end

@testset "faces_map" begin
    fmap_a, fmap_b = faces_map(a, b, Ia, Ib)
    @test fmap_a[Face(6)] == fmap_b[Face(1)]
    @test fmap_a[Face(5)] == fmap_b[Face(2)]
    @test fmap_a[Face(2)] == fmap_b[Face(3)]
    @test fmap_a[Face(3)] == fmap_b[Face(4)]
    @test fmap_a[Face(4)] == fmap_b[Face(5)]
end

@testset "contract" begin
    rzrx = contract(a, b, Ia, Ib)
    length(rzrx.half_edges_faces) == 42
end

@testset "edges" begin
    edges(a)
end

a = quon_rz()
b = quon_rx()
Ia = Vertex[9, 8, 7, 6, 5, 13]
Ib = Vertex[9, 1, 2, 3, 4, 13]

plot(a)
plot(b)

c = contract(a, b, Ia, Ib)
plot(c)