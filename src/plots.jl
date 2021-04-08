using Compose

function plot(q::QuonTait)
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
    normal_pos = Dict()
    for (k, v) in q.locations
        normal_pos[k] = ((v[1] - x_min + 0.5)/x_len, (v[2] - y_min + 0.5)/y_len)
    end

    es = []
    bds = []
    boundaries = [q.inputs; q.outputs]
    for he_id in half_edges(q)
        s = src(q, he_id)
        d = dst(q, he_id)
        if s < d
            if s in boundaries || d in boundaries
                push!(bds, [normal_pos[s], normal_pos[d]])
            else
                push!(es, [normal_pos[s], normal_pos[d]])
            end
        end
    end
    ct_bds = (context(), line(bds), stroke("green"), linewidth(1mm))
    ct_es = (context(), line(es), stroke("white"), linewidth(1mm))

    vs = vertices(q)
    vs_label = ["$v" for v in vs]
    vs_x = [normal_pos[v][1] for v in vs]
    vs_y = [normal_pos[v][2] for v in vs]
    ct_vs = (context(), text(vs_x, vs_y, vs_label), stroke("red"))

    fs = faces(q)[2:end]
    fs_label = ["$f" for f in fs]
    fs_x = []
    fs_y = []
    for f in fs
        hes_f = trace_face(q, f)
        vs_f = [src(q, he) for he in hes_f]
        push!(fs_x, sum([normal_pos[v][1] for v in vs_f])/length(vs_f))
        push!(fs_y, sum([normal_pos[v][2] for v in vs_f])/length(vs_f))
    end
    ct_fs = (context(), text(fs_x, fs_y, fs_label), stroke("blue"))

    return compose(context(), 
        ct_fs, 
        ct_vs, ct_es, ct_bds
    )
end