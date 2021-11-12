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
            if src(g, he) == dst(g, he) # self-loop
                (v1, v2) = (hes_map[he], hes_map[twin(g, he)])
                v1 < v2 && push!(new_es, (v1, v2))
            end
        else
            s = src(g, he)
            d = dst(g, he)
            s <= d && (push!(new_es, (vs_map[s], vs_map[d])))
        end
    end
    
    return v_locs, he_locs, vs_map, hes_map, new_vs, new_es
end

function plot_graph(vs, es, locs; 
        vlabel = nothing, vsize = nothing, vfill = nothing, 
        vstroke = nothing, vrots = nothing,
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
    vrots !== nothing && (vrots = [vrots[v_map_inv[v]] for v in 1:N])
    return plot_graph(adjmat, [locs[v_map[v]] for v in vs]; 
        vlabel = vlabel, vsize = vsize, vfill = vfill, vstroke = vstroke, vrots = vrots,
        tfill = tfill, tsize = tsize, kwargs...)
end

function plot_graph(adjmat, locs; padding = 0.2, 
        vlabel = nothing, vsize = nothing, vfill = nothing, 
        vstroke = nothing, vrots = nothing,
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
    vrots === nothing && (vrots = [0.0 for v in 1:N])
    line_locs = []
    for i = 1:N
        for j = i:N
            adjmat[i, j] && push!(line_locs, [(x_locs[i], y_locs[i]), (x_locs[j], y_locs[j])])
        end
    end
    
    return compose(context(units=UnitBox(x_min-padding, y_min-padding, x_max-x_min+2padding, y_max-y_min+2padding)),
        ((context(), text(x_locs[v], y_locs[v], vlabel[v], hcenter, vcenter, Rotation(vrots[v], x_locs[v], y_locs[v])), 
            fill(tfill[v]), fontsize(tsize[v])) for v in 1:N)...,
        (context(), circle(x_locs, y_locs, vsize), 
            fill(vfill), stroke(vstroke)),
        (context(), line(line_locs), stroke("gray")),
        (context(), rectangle(), fill("white"))
    )
end

function plot_planar(g::PlanarMultigraph; kwargs...)
    v_locs, he_locs, v_map, hes_map, new_vs, new_es = generate_locs(g)
    locs = Dict()
    for v in vertices(g)
        locs[v_map[v]] = v_locs[v]
    end
    for (he, v) in hes_map
        locs[v] = he_locs[he]
    end
    vlabel = Dict(v1 => "$v0" for (v0, v1) in v_map)
    return plot_graph(new_vs, new_es, locs; vlabel = vlabel, kwargs...)
end

function plot_planar(q::Tait; scale = 1, kwargs...)
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
    vrots = Dict()
    tfill = Dict()
    tsize = Dict()
    for v in vertices(q)
        locs[v_map[v]] = (v_locs[v][1] * scale, v_locs[v][2] * scale)
        vsize[v_map[v]] = is_open_vertex(q, v) ? 0.06 : 0.04
        vfill[v_map[v]] = is_open_vertex(q, v) ? "white" : (is_genus(q, v) ? "white" : "royalblue")
        vstroke[v_map[v]] = is_open_vertex(q, v) ? "white" : "royalblue"
        vrots[v_map[v]] = 0.0
        tfill[v_map[v]] = "red"
        tsize[v_map[v]] = 8/scale*pt
    end
    vlabel = Dict(v1 => "$v0" for (v0, v1) in v_map)
    for v in q.inputs
        vlabel[v_map[v]] = "$(v)ᵢ"
    end
    for v in q.outputs
        vlabel[v_map[v]] = "$(v)ₒ"
    end
    self_loops = Set{Int}()
    for (he, v) in hes_map
        if src(q, he) != dst(q, he) # not a self-loop
            locs[v] = (he_locs[he][1]*scale, he_locs[he][2]*scale)
            v_src = v_map[src(q, he)]
            v_dst = v_map[dst(q, he)]
            vrots[v] = angle(-(locs[v_src][1] - locs[v_dst][1]) + im * (locs[v_src][2] - locs[v_dst][2]))
            if is_open_half_edge(q, he)
                vfill[v] = "transparent"
                vstroke[v] = "transparent"
                vsize[v] = 0.04/scale
                vlabel[v] = "-$he→\n←$(twin(q, he))-"
            else
                vfill[v] = "gray"
                vstroke[v] = "gray"
                vsize[v] = 0.04/scale
                vlabel[v] = "-$he→\n$(quon_param(q, he))\n←$(twin(q, he))-"
            end
            tfill[v] = "black"
            tsize[v] = 5/scale*pt
        else
            he1 = he
            he2 = twin(q, he)
            if !(he1 in self_loops) 
                push!(self_loops, he1, he2)
                continue
            end
            v1 = v
            v2 = hes_map[he2]
            locs[v1] = (he_locs[he1][1]*scale, he_locs[he1][2]*scale)
            locs[v2] = (he_locs[he2][1]*scale, he_locs[he2][2]*scale)
            vfill[v1] = "gray"
            vfill[v2] = "gray"
            vstroke[v1] = "gray"
            vstroke[v2] = "gray"
            vrots[v1] = angle(-(locs[v1][1] - locs[v2][1]) + im * (locs[v1][2] - locs[v2][2]))
            vrots[v2] = 0.0
            if is_open_half_edge(q, he1)
                vfill[v1] = "transparent"
                vfill[v2] = "transparent"
                vstroke[v1] = "transparent"
                vstroke[v2] = "transparent"
                vsize[v1] = 0.04/scale
                vsize[v2] = 0.04/scale
                vlabel[v1] = "-$he→\n←$(twin(q, he))-"
                vlabel[v2] = ""
            else
                vfill[v1] = "gray"
                vfill[v2] = "gray"
                vstroke[v1] = "gray"
                vstroke[v2] = "gray"
                vsize[v1] = 0.04/scale
                vlabel[v1] = "-$he→\n$(quon_param(q, he))\n←$(twin(q, he))-"
                vsize[v2] = 0
                vlabel[v2] = ""
            end
            tfill[v1] = "black"
            tfill[v2] = "black"
            tsize[v1] = 5/scale*pt
            tsize[v2] = 5/scale*pt
        end
    end

    hes_recorded = Set(keys(hes_map))
    v_max = length(new_vs)
    for he in half_edges(q)
        (he in hes_recorded || twin(q, he) in hes_recorded) && continue
        push!(hes_recorded, he)
        v_max += 1
        push!(new_vs, v_max)
        v1 = v_map[src(q, he)]
        v2 = v_map[dst(q, he)]
        push!(new_es, (v_max, v1))
        push!(new_es, (v_max, v2))
        locs[v_max] = ((locs[v1][1] + locs[v2][1])/2, (locs[v1][2] + locs[v2][2])/2)
        tfill[v_max] = "black"
        tsize[v_max] = 5/scale*pt
        vrots[v_max] = angle(-(locs[v1][1] - locs[v2][1]) + im * (locs[v1][2] - locs[v2][2]))
        if is_open_half_edge(q, he)
            vfill[v_max] = "transparent"
            vstroke[v_max] = "transparent"
            vsize[v_max] = 0.04/scale
            vlabel[v_max] = "-$(he)→\n←$(twin(q, he))-"
        else
            vfill[v_max] = "gray"
            vstroke[v_max] = "gray"
            vsize[v_max] = 0.04/scale
            vlabel[v_max] = "-$(he)→\n$(quon_param(q, he))\n←$(twin(q, he))-"
        end
    end

    return plot_graph(new_vs, new_es, locs; 
        vlabel = vlabel, vsize = vsize, vfill = vfill, 
        vstroke = vstroke, vrots = vrots,
        tfill = tfill, tsize = tsize, kwargs...)
end