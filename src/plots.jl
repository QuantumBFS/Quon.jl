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
    for e in edges(qo.pg)
        push!(es, [normal_pos[e[1]], normal_pos[e[2]]])
    end
    ct_es = (context(), line(es), stroke("white"), linewidth(1mm))

    vs = collect(keys(qo.pos))
    vs_label = ["$(v.id)" for v in vs]
    vs_x = [normal_pos[v][1] for v in vs]
    vs_y = [normal_pos[v][2] for v in vs]
    ct_vs = (context(), text(vs_x, vs_y, vs_label), stroke("red"))
    return compose(context(), ct_vs, ct_es
    )
end