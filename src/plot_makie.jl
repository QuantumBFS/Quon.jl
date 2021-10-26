using GLMakie
using GraphMakie
using LightGraphs
using GLMakie.Colors

using GraphMakie: plot_controlpoints!, SquareGrid
g = complete_graph(3)
tangents = Dict(1 => ((1,1),(0,-1)),
                2 => ((0,1),(0,-1)),
                3 => ((0,-1),(1,0)))
tfactor = [0.5, 0.75, (0.5, 0.25)]
v_pos = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]
f, ax, p = graphplot(g; node_pos = v_pos, tangents = tangents, tfactor = tfactor,
                     edge_width = [2.0 for _ = 1:ne(g)],
                     edge_color=[colorant"red", colorant"green", colorant"blue"],
                     node_size = [10 for _ = 1:nv(g)],
                     node_color = [colorant"black" for _ = 1:nv(g)],
                     elabels="Edge ".*repr.(1:ne(g)), elabels_distance=10)
hidedecorations!(ax); hidespines!(ax); 
ax.aspect = DataAspect()
deregister_interaction!(ax, :rectanglezoom)
plot_controlpoints!(ax, p) # show control points for demonstration

function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)

function edge_hover_action(state, idx, event, axis)
    p.edge_width[][idx]= state ? 5.0 : 2.0
    p.edge_width[] = p.edge_width[] # trigger observable
end
ehover = EdgeHoverHandler(edge_hover_action)
register_interaction!(ax, :ehover, ehover)

function node_click_action(idx, args...)
    println(idx)
    p.node_color[][idx] = rand(RGB)
    p.node_color[] = p.node_color[]
end
nclick = NodeClickHandler(node_click_action)
register_interaction!(ax, :nclick, nclick)

function edge_click_action(idx, args...)
    println(idx)
    p.edge_color[][idx] = rand(RGB)
    p.edge_color[] = p.edge_color[]
end
eclick = EdgeClickHandler(edge_click_action)
register_interaction!(ax, :eclick, eclick)

function node_drag_action(state, idx, event, axis)
    p[:node_pos][][idx] = event.data
    p[:node_pos][] = p[:node_pos][]
end
ndrag = NodeDragHandler(node_drag_action)
register_interaction!(ax, :ndrag, ndrag)