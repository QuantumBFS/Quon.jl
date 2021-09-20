using PlanarMaps: PlanarMap, EmbeddedMap
using NetworkLayout: Spring, Stress, SFDP

function generate_locs(g::PlanarMultigraph; he_start = nothing)
    g_simp, hes_splitted = simple_connected_planar_graph(g)
    g_norm, v_map, he_map, f_map = normalize(g_simp)
    he0 = (he_start === nothing ? surrounding_half_edge(g_norm, f_map[0]) : he_map[he_start])
    
    adjlist = [[dst(g_norm, he) for he in trace_vertex(g_norm, v)] for v in 1:nv(g_norm)]
    p = PlanarMap(adjlist)
    e = EmbeddedMap(p; outeredge = (src(g_norm, he0), dst(g_norm, he0)))
    locs_norm = e.locs
    v_locs = Dict(v => locs_norm[v_norm] for (v, v_norm) in v_map)
    he_locs = Dict(he => locs_norm[v_map[v_he]] for (he, v_he) in hes_splitted)
    
    N = nv(g_simp)
    hes_map = Dict(he => v_map[v] for (he, v) in hes_splitted)
    vs_map = Dict(v => v_map[v] for v in vertices(g))
    new_vs = collect(1:N)
    new_es = []
    for he in half_edges(g)
        if haskey(hes_map, he)
            push!(new_es, (vs_map[src(g, he)], hes_map[he]))
        else
            s = src(g, he)
            d = dst(g, he)
            s <= d && (push!(new_es, (vs_map[s], vs_map[d])))
        end
    end
    
    return v_locs, he_locs, vs_map, hes_map, new_vs, new_es
end

function plot_graph(vs, es, locs; 
        vlabel = nothing, vsize = nothing, vfill = nothing, vstroke = nothing,
        tfill = nothing, tsize = nothing, kwargs...)
    N = length(vs)
    v_map = Dict(zip(vs, 1:N))
    v_map_inv = Dict(zip(1:N, vs))
    adjmat = zeros(Bool, N, N)
    for e in es
        adjmat[v_map[e[1]], v_map[e[2]]] = true
        adjmat[v_map[e[2]], v_map[e[1]]] = true
    end
    vlabel !== nothing && (vlabel = [vlabel[v_map_inv[v]] for v in 1:N])
    vlabel === nothing && (vlabel = ["$(v_map_inv[v])" for v in 1:N])
    vsize !== nothing && (vsize = [vsize[v_map_inv[v]] for v in 1:N])
    vfill !== nothing && (vfill = [vfill[v_map_inv[v]] for v in 1:N])
    vstroke !== nothing && (vstroke = [vstroke[v_map_inv[v]] for v in 1:N])
    tfill !== nothing && (tfill = [tfill[v_map_inv[v]] for v in 1:N])
    tsize !== nothing && (tsize = [tsize[v_map_inv[v]] for v in 1:N])
    return plot_graph(adjmat, [locs[v_map[v]] for v in vs]; 
        vlabel = vlabel, vsize = vsize, vfill = vfill, vstroke = vstroke,
        tfill = tfill, tsize = tsize, kwargs...)
end

function plot_graph(adjmat, locs; padding = 0.2, 
        vlabel = nothing, vsize = nothing, vfill = nothing, vstroke = nothing,
        tfill = nothing, tsize = nothing,
        layout = nothing)
    layout === :spring && (locs = Spring(; initialpos = locs, C = 0.2)(adjmat))
    layout === :stress && (locs = Stress(; initialpos = locs)(adjmat))
    layout === :sfdp && (locs = SFDP(; initialpos = locs, K = 0.01)(adjmat))
    
    N = size(adjmat, 1)
    x_locs = [loc[1] for loc in locs]
    y_locs = -[loc[2] for loc in locs]
    x_min = minimum(x_locs)
    x_max = maximum(x_locs)
    y_min = minimum(y_locs)
    y_max = maximum(y_locs)
    vlabel === nothing && (vlabel = ["$v" for v in 1:N])
    vsize === nothing && (vsize = [0.04 for _ in 1:N])
    vfill === nothing && (vfill = ["white" for _ in 1:N])
    vstroke === nothing && (vstroke = ["black" for _ in 1:N])
    tfill === nothing && (tfill = ["red" for _ = 1:N])
    tsize === nothing && (tsize = [8pt for _ = 1:N])
    line_locs = []
    for i = 1:N
        for j = i:N
            adjmat[i, j] && push!(line_locs, [(x_locs[i], y_locs[i]), (x_locs[j], y_locs[j])])
        end
    end
    
    return compose(context(units=UnitBox(x_min-padding, y_min-padding, x_max-x_min+2padding, y_max-y_min+2padding)),
        ((context(), text(x_locs[v], y_locs[v], vlabel[v], hcenter, vcenter), 
            fill(tfill[v]), fontsize(tsize[v])) for v in 1:N)...,
        (context(), circle(x_locs, y_locs, vsize), 
            fill(vfill), stroke(vstroke)),
        (context(), line(line_locs), stroke("gray")),
        (context(), rectangle(), fill("white"))
    )
end

function plot_planar(g::PlanarMultigraph; scale = 10)
    v_locs, he_locs, v_map, hes_map, new_vs, new_es = generate_locs(g)
    locs = Dict()
    for v in vertices(g)
        locs[v_map[v]] = v_locs[v]
    end
    for (he, v) in hes_map
        locs[v] = he_locs[he]
    end
    return plot_graph(new_vs, new_es, locs)
end

function plot_planar(q::Tait; kwargs...)
    if length(q.inputs) > 0
        hes0 = trace_vertex(q, q.inputs[1])
        he_start = hes0[findfirst(he -> (face(q, he) == 0), hes0)]
        he_start = prev(q, he_start)
    else
        he_start = nothing
    end
    v_locs, he_locs, v_map, hes_map, new_vs, new_es = generate_locs(q.g; he_start = he_start)
    locs = Dict()
    vsize = Dict()
    vfill = Dict()
    vstroke = Dict()
    tfill = Dict()
    tsize = Dict()
    for v in vertices(q)
        locs[v_map[v]] = v_locs[v]
        vsize[v_map[v]] = is_open_vertex(q, v) ? 0.06 : 0.04
        vfill[v_map[v]] = is_open_vertex(q, v) ? "white" : (is_genus(q, v) ? "white" : "royalblue")
        vstroke[v_map[v]] = is_open_vertex(q, v) ? "white" : "royalblue"
        tfill[v_map[v]] = "red"
        tsize[v_map[v]] = 8pt
    end
    vlabel = Dict(v1 => "$v0" for (v0, v1) in v_map)
    for v in q.inputs
        vlabel[v_map[v]] = "$(v)ᵢ"
    end
    for v in q.outputs
        vlabel[v_map[v]] = "$(v)ₒ"
    end
    for (he, v) in hes_map
        locs[v] = he_locs[he]
        vsize[v] = 0.04
        vfill[v] = "gray"
        vstroke[v] = "gray"
        vlabel[v] = "$(phase(q, he))"
        tfill[v] = "black"
        tsize[v] = 5pt
    end

    hes_recorded = Set(keys(hes_map))
    v_max = length(new_vs)
    for he in half_edges(q)
        (he in hes_recorded || twin(q, he) in hes_recorded) && continue
        
        v_max += 1
        push!(new_vs, v_max)
        v1 = v_map[src(q, he)]
        v2 = v_map[dst(q, he)]
        push!(new_es, (v_max, v1))
        push!(new_es, (v_max, v2))
        locs[v_max] = ((locs[v1][1] + locs[v2][1])/2, (locs[v1][2] + locs[v2][2])/2)
        tfill[v_max] = "black"
        tsize[v_max] = 5pt
        vfill[v_max] = "gray"
        vstroke[v_max] = "gray"
        if is_open_half_edge(q, he)
            vsize[v_max] = 0.0
            vlabel[v_max] = ""
        else
            vsize[v_max] = 0.04
            vlabel[v_max] = "$(phase(q, he))"
        end
    end

    return plot_graph(new_vs, new_es, locs; 
        vlabel = vlabel, vsize = vsize, vfill = vfill, vstroke = vstroke,
        tfill = tfill, tsize = tsize, kwargs...)
end