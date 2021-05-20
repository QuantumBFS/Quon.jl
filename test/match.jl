using Test
using Quon

@testset "string genus" begin
    c = contract!(tait_copy(), tait_copy(), [2, 3], [1, 3])
    m = Quon.match(Quon.Rule{:string_genus}(), c)
    @test length(m) == 1
    @test m[1].vertices[1] == 4
end
