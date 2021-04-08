using Quon
using Test

@testset "yang_baxter.jl" begin
    include("utils/yang_baxter.jl")
end

@testset "graphs/tait.jl" begin
    include("graphs/tait.jl")
end

@testset "plots.jl" begin
    include("plots.jl")
end