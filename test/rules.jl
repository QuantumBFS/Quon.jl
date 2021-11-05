using Test
using Quon
using Quon: check

@testset "string-genus Rule" begin
    c = contract!(tait_copy(), tait_copy(), [2, 3], [1, 3])
    m = Quon.match(Quon.Rule{:string_genus}(), c)
    @test length(m) == 1
    @test check(c, m[1])
    g = m[1].vertices[1]
    rewrite!(c, m[1])
    @test !(g in c.genuses)
end

@testset "perm_rz" begin
    rz2 = contract!(tait_rz(0*im), tait_rz(0*im))
    m = Quon.match(Quon.Rule{:perm_rz}(), rz2)
    @test length(m) == 2
    @test check(rz2, m[1]) && check(rz2, m[2])

    rz = tait_rz(0*im)
    m = Quon.match(Quon.Rule{:perm_rz}(), rz)
    @test check(rz, m[1])
    Quon.rewrite!(rz, m[1])
    @test rz.g.f_max == 4
end

@testset "Yang-Baxter Rule" begin 
    ry = contract!(contract!(tait_rz(pi/2*im), tait_rx(-π/2*im)), tait_rz(π/2*im))
    m1 = match(Rule{:yang_baxter_triangle}(), ry)
    @test length(m1) == 1
    @test check(ry, m1[1])
    rewrite!(ry, m1[1])

    m2 = match(Rule{:yang_baxter_star}(), ry)
    @test length(m2) == 1
    @test check(ry, m2[1])
    Quon.rewrite!(ry, m2[1])
    m3 = match(Rule{:yang_baxter_triangle}(), ry)
    @test length(m3) == 1
end

@testset "Fusion Rule" begin
    rx2 = contract!(tait_rx(π/2*im), tait_rx(π/2*im))
    mx = match(Rule{:x_fusion}(), rx2)
    @test check(rx2, mx[1])
    rewrite!(rx2, mx[1])
    @test length(rx2.quon_params) == 2

    rz2 = contract!(tait_rz(π/2*im), tait_rz(π/2*im))
    mz = match(Rule{:z_fusion}(), rz2)
    @test check(rz2, mz[1])
    rewrite!(rz2, mz[1])
    @test length(rz2.quon_params) == 2
end

@testset "charge remove" begin
    rx2 = contract!(tait_rx(pi*im), tait_rx(pi*im))
    m_charge_rm_v = match(Rule(:charge_rm_v), rx2)
    @test check(rx2, m_charge_rm_v[1])
    rewrite!(rx2, m_charge_rm_v[1])
    @test nv(rx2) == 5

    rz2 = contract!(tait_rz(pi*im), tait_rz(pi*im))
    m_charge_rm_f = match(Rule(:charge_rm_f), rz2)
    @test check(rz2, m_charge_rm_f[1])
    rewrite!(rz2, m_charge_rm_f[1])
    @test nv(rz2) == 5
end

@testset "SWAP-genus Rule" begin
    s2 = contract!(tait_swap(), tait_swap())
    m_swap_genus = match(Rule(:swap_genus), s2)
    @test check(s2, m_swap_genus[1])
    rewrite!(s2, m_swap_genus[1])
    @test !check(s2, m_swap_genus[1])
end

# TODO: fix locations for cz
# @testset "CZ^2 = I" begin
#     cz = tait_cz()
#     contract!(cz, tait_cz())
#     m = match(Rule{:perm_rz}(), cz)
#     Quon.rewrite!(cz, m[3])
#     Quon.rewrite!(cz, m[9])
#     m_string_genus = match(Quon.Rule{:string_genus}(), cz)
#     @test length(m_string_genus) == 1
#     Quon.rewrite!(cz, m_string_genus[1])
#     m_fusion = match(Quon.Rule{:z_fusion}(), cz)
#     @test length(m_fusion) == 3
#     Quon.rewrite!(cz, m_fusion[1])
#     Quon.rewrite!(cz, m_fusion[2])
#     Quon.rewrite!(cz, m_fusion[3])
#     m_yb = match(Quon.Rule{:yang_baxter_triangle}(), cz)
#     @test length(m_yb) == 1
#     rewrite!(cz, m_yb[1])
#     m_id = match(Quon.Rule{:identity}(), cz)
#     Quon.rewrite!(cz, m_id[1])
#     Quon.rewrite!(cz, m_id[3])
#     Quon.rewrite!(cz, m_id[5])
#     @test length(m_id) == 6
#     m_gf = match(Rule(:genus_fusion), cz)
#     @test length(m_gf) == 2
#     rewrite!(cz, m_gf[2])
# end