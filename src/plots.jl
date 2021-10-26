using YaoPlots
include("plot_compose.jl")
include("planar_maps.jl")

function YaoPlots.plot(q::Tait; backend = :planar, kwargs...)
    backend === :compose && return plot_compose(q; kwargs...)
    backend === :planar && return plot_planar(q; kwargs...)
end

function YaoPlots.plot(g::PlanarMultigraph; backend = :planar, kwargs...)
    backend === :planar && return plot_planar(g; kwargs...)
end