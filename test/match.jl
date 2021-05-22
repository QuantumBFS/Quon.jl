using Test
using Quon

@testset "string genus" begin
    c = contract!(tait_copy(), tait_copy(), [2, 3], [1, 3])
    m = Quon.match(Quon.Rule{:string_genus}(), c)
    @test length(m) == 1
    @test m[1].vertices[1] == 4
end

@testset "perm_rz" begin
    rz2 = contract!(tait_rz(0*im), tait_rz(0*im))
    m = Quon.match(Quon.Rule{:perm_rz}(), rz2)
    @test length(m) == 2

    rz = tait_rz(0*im)
    m = Quon.match(Quon.Rule{:perm_rz}(), rz)
    Quon.rewrite!(rz, m[1])
    @test rz.g.f_max == 4
end

@testset "Yang-Baxter Rule" begin 
    ry = contract!(contract!(tait_rz(0.0*im), tait_rx(-π/2*im)), tait_rz(π/2*im))
    m1 = match(Rule{:yang_baxter_triangle}(), ry)
    @test length(m1) == 1
    rewrite!(ry, m1[1])

    m2 = match(Rule{:yang_baxter_star}(), ry)
    @test length(m2) == 1
    Quon.rewrite!(ry, m2[1])
    m3 = match(Rule{:yang_baxter_triangle}(), ry)
    @test length(m3) == 1
end