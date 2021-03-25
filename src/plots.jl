export plot
using Compose

function plot(qo::QuonGraph)
    x_min, x_max = (0.0, 1.0)
    y_min, y_max = (0.0, 1.0)
    for v in values(qo.pos)
        v[1] < x_min && (x_min = v[1])
        v[1] > x_max && (x_max = v[1])
        v[2] < y_min && (y_min = v[2])
        v[2] > y_max && (y_max = v[2])
    end
    x_len = x_max - x_min + 1
    y_len = y_max - y_min + 1
    normal_pos = Dict()
    for (k, v) in qo.pos
        normal_pos[k] = ((v[1] - x_min + 0.5)/x_len, (v[2] - y_min + 0.5)/y_len)
    end

    es = []
    bds = []
    for he in half_edges(qo.pg)
        if src(he).id < dst(he).id 
            if !is_boundary(qo.pg, he) && !is_boundary(qo.pg, twin(he))
                push!(es, [normal_pos[src(he)], normal_pos[dst(he)]])
            else
                push!(bds, [normal_pos[src(he)], normal_pos[dst(he)]])
            end
        end
    end
    ct_bds = (context(), line(bds), stroke("green"), linewidth(1mm))
    ct_es = (context(), line(es), stroke("white"), linewidth(1mm))

    vs = collect(keys(qo.pos))
    vs_label = ["$(v.id)" for v in vs]
    vs_x = [normal_pos[v][1] for v in vs]
    vs_y = [normal_pos[v][2] for v in vs]
    ct_vs = (context(), text(vs_x, vs_y, vs_label), stroke("red"))

    fs = faces(qo.pg)
    deleteat!(fs, findfirst(isequal(Face(0)), fs))
    fs_label = ["$(f.id)" for f in fs]
    fs_x = []
    fs_y = []
    for f in fs
        vs_f = trace_face(qo.pg, f)[2]
        push!(fs_x, sum([normal_pos[v][1] for v in vs_f])/length(vs_f))
        push!(fs_y, sum([normal_pos[v][2] for v in vs_f])/length(vs_f))
    end
    ct_fs = (context(), text(fs_x, fs_y, fs_label), stroke("blue"))

    return compose(context(), ct_fs, ct_vs, ct_es, ct_bds
    )
end