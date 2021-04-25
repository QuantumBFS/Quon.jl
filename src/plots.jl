using Compose

function plot(q::Tait; show_half_edges = true, show_faces = true, face_color = "salmon", background = "white")
    x_min, x_max = (0.0, 1.0)
    y_min, y_max = (0.0, 1.0)
    for v in values(q.locations)
        v[1] < x_min && (x_min = v[1])
        v[1] > x_max && (x_max = v[1])
        v[2] < y_min && (y_min = v[2])
        v[2] > y_max && (y_max = v[2])
    end
    x_len = x_max - x_min + 1
    y_len = y_max - y_min + 1
    original_pos = q.locations
    normal_pos = Dict()
    for (k, v) in q.locations
        normal_pos[k] = ((v[1] - x_min + 0.5)/x_len, (v[2] - y_min + 0.5)/y_len)
    end

    es = []
    bds = []
    boundaries = [q.inputs; q.outputs]
    hes_label = String[]
    hes_x = Float64[]
    hes_y = Float64[]
    hes_rot = Float64[]
    for he_id in half_edges(q)
        s = src(q, he_id)
        d = dst(q, he_id)
        pos_s = normal_pos[s]
        pos_d = normal_pos[d]
        if s < d
            if s in boundaries || d in boundaries
                push!(bds, [pos_s, pos_d])
            else
                push!(es, [pos_s, pos_d])
            end
        end
        push!(hes_label, "$he_id\n→")
        push!(hes_x, (pos_s[1]+pos_d[1])/2)
        push!(hes_y, (pos_s[2]+pos_d[2])/2)
        push!(hes_rot, angle((original_pos[d][1]-original_pos[s][1])+im*(original_pos[d][2]-original_pos[s][2])))
    end
    
    ct_bds = (context(), line(bds), stroke("black"), linewidth(1mm))
    ct_es = (context(),
        (context(), line(es), stroke("white"), linewidth(0.8mm)),
        (context(), line(es), stroke("black"), linewidth(1.2mm)),
    )
    rot = Rotation.(hes_rot, hes_x, hes_y)
    ct_hes = (context(), 
        text(hes_x, hes_y, hes_label, [hcenter], [vbottom], rot), 
        fill("gray"),
        fontsize(8pt)
    )

    vs = vertices(q)
    vs_label = ["$v" for v in vs]
    vs_x = [normal_pos[v][1] for v in vs]
    vs_y = [normal_pos[v][2] for v in vs]
    ct_vs = context()
    for i in 1:length(vs)
        ct_v = (context(vs_x[i]*w-0.5cm, vs_y[i]*h-0.5cm, 1cm, 1cm),
            is_genus(q, vs[i]) ? (context(), arc([0.5], [0.4, 0.8], [0.4cm, 0.2cm], [0.5, pi+0.5], [pi-0.5, -0.5], [false]), fill("transparent"), stroke("black"), linewidth(0.5mm)) : context(),
            (context(), text(0.5, 0.5, vs_label[i], hcenter, is_genus(q, vs[i]) ? vbottom : vcenter), fill("red"), font("Helvetica-Bold")),

            (context(), circle(), is_open(q, vs[i]) ? fill("black") : fill("white"), stroke("black"), linewidth(0.8mm)),
        )
        ct_vs = compose(context(), ct_vs, ct_v)
    end

    fs = faces(q)[2:end]
    fs_label = ["$f" for f in fs]
    fs_x = []
    fs_y = []
    fs_points = []
    for f in fs
        hes_f = trace_face(q, f)
        vs_f = [src(q, he) for he in hes_f]
        push!(fs_x, sum([normal_pos[v][1] for v in vs_f])/length(vs_f))
        push!(fs_y, sum([normal_pos[v][2] for v in vs_f])/length(vs_f))
        push!(fs_points, [(normal_pos[v][1], normal_pos[v][2]) for v in vs_f])
    end
    ct_fs = (context(), text(fs_x, fs_y, fs_label, [hcenter], [vcenter]), fill("blue"), font("Helvetica-Bold"))
    ct_polys = context()
    for ps in fs_points
        ct_polys = compose(context(), 
            ct_polys,
            (context(), polygon(ps), fill(face_color), 
                # fillopacity(0.3)
            )
        )
    end

    bg = (context(), rectangle(), fill(background))
    set_default_graphic_size(3x_len*cm, 3y_len*cm)
    return compose(context(), 
        show_faces ? ct_fs : context(), 
        show_half_edges ? ct_hes : context(),
        ct_vs, 
        ct_es, 
        ct_bds,
        ct_polys,
        bg
    )
end