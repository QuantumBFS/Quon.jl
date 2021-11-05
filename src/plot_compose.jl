using Compose

function plot_compose(q::Tait; show_string = true, show_tait = true, show_half_edges = true, show_faces = true, face_color = "salmon", background = "white", 
        start_hes = Dict(), radius = 0.5)
    x_len, y_len, normal_pos = generate_normalized_locations(q)
    es, bds, hes_label, hes_x, hes_y, hes_rot = generate_half_edges_info(q, normal_pos)
    ctrl_points = generate_ctrl_points(q, normal_pos, radius, start_hes)
    update_half_edge_info!(q, hes_x, hes_y, hes_rot, normal_pos, ctrl_points)
    ct_edges = generate_context_edges(q, es, bds, hes_label, hes_x, hes_y, hes_rot, ctrl_points, show_half_edges)
    ct_vs, ct_genus = generate_context_vertices(q, normal_pos, show_tait)
    ct_fs = generate_context_faces(q, normal_pos, ctrl_points, show_faces)
    
    bg = (context(), rectangle(), fill(background))
    set_default_graphic_size(3(x_len+1)*cm, 3(y_len+1)*cm)
    ct_string = generate_string(q, hes_x, hes_y, hes_rot, radius/4)
    return compose(
        context(units=UnitBox(0, 0, x_len, y_len; 
            leftpad=1.5cm, rightpad=1.5cm, toppad=1.5cm, bottompad=1.5cm)), 
        show_string ? ct_string : context(),
        ct_genus,
        show_tait ? ct_vs : context(),
        show_tait ? ct_edges : context(), 
        show_tait ? ct_fs : context(),
        bg
    )
end

# function generate_locations(q::Tait)
#     vs_boundary = [q.inputs; q.outputs; genuses(q)]
#     vs_proceeded
    
#     frontier = copy(q.inputs)
#     while !isempty(frontier)

#     end
# end

function generate_normalized_locations(q::Tait)
    x_min, x_max = (0.0, 1.0)
    y_min, y_max = (0.0, 1.0)
    for v in values(q.locations)
        v[1] < x_min && (x_min = v[1])
        v[1] > x_max && (x_max = v[1])
        v[2] < y_min && (y_min = v[2])
        v[2] > y_max && (y_max = v[2])
    end
    x_len = x_max - x_min
    y_len = y_max - y_min
    original_pos = q.locations
    normal_pos = Dict()
    for (k, v) in q.locations
        normal_pos[k] = ((v[1] - x_min), (v[2] - y_min))
    end
    return x_len, y_len, normal_pos
end

function generate_half_edges_info(q::Tait, normal_pos)
    es = Dict()
    bds = Dict()
    boundaries = [q.inputs; q.outputs]
    hes_label = Dict()
    hes_x = Dict()
    hes_y = Dict()
    hes_rot = Dict()

    for he_id in half_edges(q)
        s = src(q, he_id)
        d = dst(q, he_id)
        pos_s = normal_pos[s]
        pos_d = normal_pos[d]
        if s in boundaries || d in boundaries
            bds[he_id] = pos_s
        else
            es[he_id] = pos_s
        end
        hes_label[he_id] = "$he_id\n→"
        hes_x[he_id] = (pos_s[1]+pos_d[1])/2
        hes_y[he_id] = (pos_s[2]+pos_d[2])/2
        hes_rot[he_id] = angle((normal_pos[d][1]-normal_pos[s][1])+im*(normal_pos[d][2]-normal_pos[s][2]))
    end
    return es, bds, hes_label, hes_x, hes_y, hes_rot
end

function generate_ctrl_points(q::Tait, normal_pos, radius, start_hes)
    ctrl_points = Dict()
    for v in vertices(q)
        if !haskey(start_hes, v)
            hes = trace_vertex(q, v)
            n = length(hes)
            θ = -π/2
            pos_v = normal_pos[v]
            for i in 1:n
                he = hes[i]
                θ += 2π/n
                ctrl_points[he] = (pos_v[1] + cos(θ)*radius, pos_v[2] - sin(θ)*radius)
            end
        else
            start_he, θ = start_hes[v]
            hes_v = trace_vertex(q, v)
            n = length(hes_v)
            thetas = [2π*i/n for i = 1:n]
            start_id = findfirst(isequal(start_he), hes_v)
            thetas .+= θ - thetas[start_id]
            pos_v = normal_pos[v]
            for i = 1:n
                he = hes_v[i]
                ctrl_points[he] = (pos_v[1] + cos(thetas[i])*radius, pos_v[2] - sin(thetas[i])*radius)
            end
        end
    end
    return ctrl_points
end

function update_half_edge_info!(q::Tait, hes_x, hes_y, hes_rot, normal_pos, ctrl_points)
    for he in half_edges(q)
        hes_x[he] = hes_x[he]/4 + 3/8*(ctrl_points[he][1] + ctrl_points[twin(q, he)][1])
        hes_y[he] = hes_y[he]/4 + 3/8*(ctrl_points[he][2] + ctrl_points[twin(q, he)][2])
        s = src(q, he)
        d = dst(q, he)
        
        hes_rot[he] = angle(
            (normal_pos[d][1] + ctrl_points[twin(q, he)][1] - normal_pos[s][1] - ctrl_points[he][1])
                + im * (normal_pos[d][2] + ctrl_points[twin(q, he)][2] - normal_pos[s][2] - ctrl_points[he][2])
        )
    end
    return 
end

function generate_context_edges(q::Tait, es, bds, hes_label, hes_x, hes_y, hes_rot, ctrl_points, show_half_edges)
    bd_hes = Set{Int}()
    for he in keys(bds)
        !(twin(q, he) in bd_hes) && push!(bd_hes, he)
    end
    bd_hes = sort(collect(bd_hes))
    inner_hes = Set{Int}()
    for he in keys(es)
        !(twin(q, he) in inner_hes) && push!(inner_hes, he)
    end
    inner_hes = sort(collect(inner_hes))
    ct_bds = length(bd_hes) > 0 ? (context(), 
        curve([bds[he] for he in bd_hes], 
            [ctrl_points[he] for he in bd_hes], 
            [ctrl_points[twin(q, he)] for he in bd_hes],
            [bds[twin(q, he)] for he in bd_hes],
        ), 
        stroke("gray"), linewidth(1mm)
    ) : context()
    ct_es = !isempty(inner_hes) ? (context(),
        (context(), 
            curve([es[he] for he in inner_hes], 
                [ctrl_points[he] for he in inner_hes], 
                [ctrl_points[twin(q, he)] for he in inner_hes],
                [es[twin(q, he)] for he in inner_hes],
            ), stroke("white"), linewidth(0.6mm)),
        (context(), 
            curve([es[he] for he in inner_hes], 
                [ctrl_points[he] for he in inner_hes], 
                [ctrl_points[twin(q, he)] for he in inner_hes],
                [es[twin(q, he)] for he in inner_hes],
            ), stroke("gray"), linewidth(1.2mm)),
    ) : context()
    rot = [Rotation(hes_rot[he], hes_x[he], hes_y[he]) for he in half_edges(q)]
    ct_hes = (context(), 
        text([hes_x[he] for he in half_edges(q)], [hes_y[he] for he in half_edges(q)], 
            [hes_label[he] for he in half_edges(q)], [hcenter], [vbottom], rot), 
        fill("gray"),
        fontsize(8pt)
    )
    ct_edges = compose(context(), 
        show_half_edges ? ct_hes : context(),
        ct_es, 
        ct_bds
    )
    return ct_edges
end

function generate_context_vertices(q, normal_pos, show_tait)
    vs = vertices(q)
    vs_label = ["$v" for v in vs]
    vs_x = [normal_pos[v][1] for v in vs]
    vs_y = [normal_pos[v][2] for v in vs]
    ct_vs = context()
    ct_genus = context()
    for i in 1:length(vs)
        if is_isolated(q, vs[i])
            ct_v = compose(
                context(vs_x[i]-0.1, vs_y[i]-0.1, 0.2, 0.2),
                (
                    context(units = UnitBox(2,2)),
                    is_genus(q, vs[i]) ? (context(), arc([1], [0.8, 1.6], [0.8, 0.4], [0.5, pi+0.5], [pi-0.5, -0.5], [false]), fill("transparent"), stroke("royalblue"), linewidth(0.5mm)) : context(),
                    (context(), text(1, 1, vs_label[i], hcenter, is_genus(q, vs[i]) ? vbottom : vcenter), fill("red"), font("Helvetica-Bold")),
                    (context(), circle(1, 1, 1), is_open_vertex(q, vs[i]) ? fill("royalblue") : fill("white"), 
                        (is_genus(q, vs[i]) && !show_tait) ? stroke("black") : stroke("royalblue"), linewidth(1)
                    )
                )
            )
        else
            ct_v = compose(
                context(vs_x[i]-0.15, vs_y[i]-0.15, 0.3, 0.3),
                (
                    context(units = UnitBox(2,2)),
                    is_genus(q, vs[i]) ? (context(), arc([1], [0.8, 1.6], [0.8, 0.4], [0.5, pi+0.5], [pi-0.5, -0.5], [false]), fill("transparent"), stroke("royalblue"), linewidth(0.5mm)) : context(),
                    (context(), text(1, 1, vs_label[i], hcenter, is_genus(q, vs[i]) ? vbottom : vcenter), fill("red"), font("Helvetica-Bold")),
                    (context(), circle(1, 1, 1), is_open_vertex(q, vs[i]) ? fill("royalblue") : fill("white"), 
                        (is_genus(q, vs[i]) && !show_tait) ? stroke("white") : stroke("royalblue"), linewidth(0.5)
                    )
                )
            )
        end
        if is_genus(q, vs[i])
            ct_genus = compose(context(), ct_genus, ct_v)
        else
            ct_vs = compose(context(), ct_vs, ct_v)
        end
    end
    return ct_vs, ct_genus
end

function generate_context_faces(q::Tait, normal_pos, ctrl_points, show_faces)
    fs = faces(q)[2:end]
    fs_label = ["$f" for f in fs]
    fs_x = []
    fs_y = []
    ct_bezi = compose(context())
    for f in fs
        hes_f = trace_face(q, f)
        vs_f = [src(q, he) for he in hes_f]
        push!(fs_x, (sum([(ctrl_points[he][1] + ctrl_points[twin(q, he)][1])/2 for he in hes_f])*0.5 + sum([normal_pos[v][1] for v in vs_f])*0.5)/(length(vs_f)))
        push!(fs_y, (sum([(ctrl_points[he][2] + ctrl_points[twin(q, he)][2])/2 for he in hes_f])*0.5 + sum([normal_pos[v][2] for v in vs_f])*0.5)/(length(vs_f)))
        start_point = normal_pos[src(q, hes_f[1])]
        bz_ctrls = [[ctrl_points[he], ctrl_points[twin(q, he)], normal_pos[dst(q, he)]] for he in hes_f]
        ct_bezi_f = compose(context(), 
            bezigon(start_point, bz_ctrls), fill("salmon"), fillopacity(0.3),
        )
        ct_bezi = compose(context(), ct_bezi, ct_bezi_f)
    end
    if length(fs_x) > 0
        ct_fs = (context(), text(fs_x, fs_y, fs_label, [hcenter], [vcenter]), fill("blue"), font("Helvetica-Bold"))
    else
        ct_fs = context()
    end
    return compose(context(),
        show_faces ? ct_fs : context(),
        show_faces ? ct_bezi : context(),
    )
end

function generate_string(q::Tait, hes_x, hes_y, hes_rot, radius)
    ct_string = compose(context())
    ct_phase = compose(context())
    for f in faces(q)
        f == 0 && continue
        he_start = surrounding_half_edge(q, f)
        he_curr = he_start
        he_next = next(q, he_curr)
        ct_radius = 3radius
        while true
            if !(is_open_vertex(q, dst(q, he_curr)) && is_open_vertex(q, src(q, he_next)))
                if is_open_vertex(q, src(q, he_curr))
                    pos_curr = (hes_x[he_curr] + radius*cos(hes_rot[he_curr]-π/2), hes_y[he_curr] + radius*sin(hes_rot[he_curr]-π/2))
                    ctrl_curr = (pos_curr[1] + ct_radius*cos(hes_rot[he_curr]), pos_curr[2] + ct_radius*sin(hes_rot[he_curr]))
                else
                    pos_curr = (hes_x[he_curr], hes_y[he_curr])
                    ctrl_curr = (pos_curr[1] + ct_radius*cos(hes_rot[he_curr] - π/4), pos_curr[2] + ct_radius*sin(hes_rot[he_curr] - π/4))
                end
                if is_open_vertex(q, dst(q, he_next))
                    pos_next = (hes_x[he_next] + radius*cos(hes_rot[he_next]-π/2), hes_y[he_next] + radius*sin(hes_rot[he_next]-π/2))
                    ctrl_next = (pos_next[1] + ct_radius*cos(hes_rot[he_next]-π), pos_next[2] + ct_radius*sin(hes_rot[he_next]-π))
                else
                    pos_next = (hes_x[he_next], hes_y[he_next])
                    ctrl_next = (pos_next[1] + ct_radius*cos(hes_rot[he_next] - 3π/4), pos_next[2] + ct_radius*sin(hes_rot[he_next] - 3π/4))
                end
                ct_string = compose(
                    ct_string,
                    (context(), 
                    curve(
                        [pos_curr], [ctrl_curr], 
                        [ctrl_next], [pos_next],
                    ),
                    stroke("black"), linewidth(1mm))
                )
            end
            he_curr = he_next
            he_next = next(q, he_curr)
            (he_curr == he_start) && break
        end
    end
    hes = Set(half_edges(q))
    while !isempty(hes)
        he = pop!(hes)
        if !is_open_vertex(q, src(q, he)) && !is_open_vertex(q, dst(q, he))
            θ = hes_rot[he]
            crx = generate_phase(hes_x[he], hes_y[he], θ, quon_param(q, he))
            ct_phase = compose(context(), ct_phase, crx)
        end
        delete!(hes, twin(q, he))
    end
    return compose(context(), ct_phase, ct_string)
end

function generate_phase(x, y, θ, p)
    if (isapprox(p.param, π/2*im; atol = quon_atol) && !is_parallel(p)) ||
        (isapprox(p.param, -π/2*im; atol = quon_atol) && is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), line([(cos(θ+π/4), sin(θ+π/4)), (cos(θ+5π/4), sin(θ+5π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))),
        )
    elseif (isapprox(p.param, -π/2*im; atol = quon_atol) && !is_parallel(p)) ||
        (isapprox(p.param, π/2*im; atol = quon_atol) && is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), line([(cos(θ+3π/4), sin(θ+3π/4)), (cos(θ+7π/4), sin(θ+7π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))),
        )
    elseif (isapprox(real(p.param), 0; atol = quon_atol) && 
        isapprox(rem(imag(p.param), 2π, RoundDown), π; atol = quon_atol) && is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), circle(cos(θ+π/2)/sqrt(2), sin(θ+π/2)/sqrt(2), 0.4), fill("black")),
                (context(), circle(cos(θ+3π/2)/sqrt(2), sin(θ+3π/2)/sqrt(2), 0.4), fill("black")),
                (context(), line([(cos(θ+π/4), sin(θ+π/4)), (cos(θ+3π/4), sin(θ+3π/4))]), stroke("black"), linewidth(1)),
                (context(), line([(cos(θ+7π/4), sin(θ+7π/4)), (cos(θ+5π/4), sin(θ+5π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))
            ),
        )
    elseif (isapprox(real(p.param), 0; atol = quon_atol) && 
            isapprox(rem(imag(p.param), 2π, RoundDown), π; atol = quon_atol) && !is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), circle(cos(θ)/sqrt(2), sin(θ)/sqrt(2), 0.4), fill("black")),
                (context(), circle(cos(θ+π)/sqrt(2), sin(θ+π)/sqrt(2), 0.4), fill("black")),
                (context(), line([(cos(θ+3π/4), sin(θ+3π/4)), (cos(θ+5π/4), sin(θ+5π/4))]), stroke("black"), linewidth(1)),
                (context(), line([(cos(θ+7π/4), sin(θ+7π/4)), (cos(θ+π/4), sin(θ+π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))
            ),
        )
    elseif (isapprox(real(p.param), 0; atol = quon_atol) && 
        isapprox(rem(imag(p.param), 2π, RoundDown), 0; atol = quon_atol) && is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), line([(cos(θ+π/4), sin(θ+π/4)), (cos(θ+3π/4), sin(θ+3π/4))]), stroke("black"), linewidth(1)),
                (context(), line([(cos(θ+7π/4), sin(θ+7π/4)), (cos(θ+5π/4), sin(θ+5π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))
            ),
        )
    elseif (isapprox(real(p.param), 0; atol = quon_atol) && 
            isapprox(rem(imag(p.param), 2π, RoundDown), 0; atol = quon_atol) && !is_parallel(p))
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), line([(cos(θ+3π/4), sin(θ+3π/4)), (cos(θ+5π/4), sin(θ+5π/4))]), stroke("black"), linewidth(1)),
                (context(), line([(cos(θ+7π/4), sin(θ+7π/4)), (cos(θ+π/4), sin(θ+π/4))]), stroke("black"), linewidth(1)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))
            ),
        )
    else
        return compose(
            context(x-0.1, y-0.1, 0.2, 0.2),
            (context(units=UnitBox(-1, -1, 2, 2)),
                (context(), text(0, 0, "$p", hcenter, vcenter, Rotation(θ, 0, 0)), fontsize(1.5)),
                (context(), circle(0, 0, 1), fill("white"), stroke("black"))
            ),
        )
    end
end