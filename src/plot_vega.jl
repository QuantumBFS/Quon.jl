using Vega
using DataFrames

function generate_d_vs(q::Tait, x_min, y_min)
    ids = sort!(collect(vertices(q)))
    vts = [Quon.is_open_vertex(q, v) ? "B" : (Quon.is_genus(q, v) ? "G" : "N") for v in ids]
    xs = [(q.locations[v][1]-x_min)*100 for v in ids]
    ys = [(q.locations[v][2]-y_min)*100 for v in ids]
    degs = [length(trace_vertex(q, v)) for v in ids]
    shifts = zeros(length(ids))
    d_vs = DataFrame(id = ids, v_type = vts, x = xs, y = ys, shift = shifts, degree = degs)
    return d_vs
end

function generate_d_edges(q::Tait)
    hes = sort!(collect(half_edges(q)))
    fs = [face(q, he) for he in hes]
    srcs = [src(q, he) for he in hes]
    dsts = [dst(q, he) for he in hes]
    rs1 = Dict{Int, Float64}()
    rs2 = Dict{Int, Float64}()
    for v in vertices(q)
        hes_v = trace_vertex(q, v)
        n = length(hes_v)
        for i in eachindex(hes_v)
            he = hes_v[i]
            rs1[he] = -(i-1)*2pi/n
            rs2[twin(q, he)] = -(i-1)*2pi/n
        end
    end
    rot_srcs = [rs1[he] for he in hes]
    rot_dsts = [rs2[he] for he in hes]
    d_edges = DataFrame(id = hes, src = srcs, dst = dsts, f = fs, rot_src = rot_srcs, rot_dst = rot_dsts)
    return d_edges
end

function generate_frame(q::Tait)
    x_min, x_max = (0.0, 1.0)
    y_min, y_max = (0.0, 1.0)
    for v in vertices(q)
        loc = q.locations[v]
        loc[1] > x_max && (x_max = loc[1])
        loc[1] < x_min && (x_min = loc[1])
        loc[2] > y_max && (y_max = loc[2])
        loc[2] < y_min && (y_min = loc[2])
    end
    return x_min, x_max, y_min, y_max
end

function plot_vega(q::Tait)
    x_min, x_max, y_min, y_max = generate_frame(q)
    h, w = (x_max-x_min)*100, (y_max-y_min)*100
    println(h, w)
    d_edges = generate_d_edges(q)
    d_vs = generate_d_vs(q, x_min, y_min)
    return @vgplot(
        $schema="https://vega.github.io/schema/vega/v5.json",
        height=h,
        width=w,
        padding=50,
        autosize="none",
        data=[
            {
                name="vertices",
                values = d_vs,
                on=[
                    {
                        modify="whichSymbol",
                        values="newLoc && {x: newLoc.x, y: newLoc.y}",
                        trigger="newLoc"
                    },
                    {
                        modify="whichSymbol",
                        values="Clicked && {shift: whichSymbol.shift + 3.14/whichSymbol.degree}",
                        trigger="Clicked"
                    }
                ],
                transform=[
                    {
                        as="shape",
                        expr="'circle'",
                        type="formula"
                    },
                    {
                        as="color",
                        expr="datum.v_type === 'B' ? '#389826' : (datum.v_type === 'N' ? '#CB3C33' : (datum.v_type === 'G' ? 'white' : 'white'))",
                        type="formula"
                    }
                ]
            },
            {
                name="edges",
                values=d_edges,
                transform=[
                    {
                        key="id",
                        fields=[
                            "dst",
                            "src"
                        ],
                        as=[
                            "t",
                            "s"
                        ],
                        from="vertices",
                        type="lookup"
                    },
                    {
                        as="c1",
                        expr="{'x': datum.s.x + ctRadius * cos(datum.rot_src + datum.s.shift), 'y': datum.s.y + ctRadius * sin(datum.rot_src + datum.s.shift)}",
                        type="formula"
                    },
                    {
                        as="c2",
                        expr="{'x': datum.t.x + ctRadius * cos(datum.rot_dst + datum.t.shift), 'y': datum.t.y + ctRadius * sin(datum.rot_dst + datum.t.shift)}",
                        type="formula"
                    },
                    {
                        as="m",
                        expr="{'x': datum.s.x*1/8 + datum.t.x*1/8 + datum.c1.x*3/8 + datum.c2.x*3/8, 'y': datum.s.y*1/8 + datum.t.y*1/8 + datum.c1.y*3/8 + datum.c2.y*3/8}",
                        type="formula"
                    },
                    {
                        as="mypath",
                        expr="join(['M', toString(datum.s.x), toString(datum.s.y), 'C', toString(datum.c1.x), toString(datum.c1.y), toString(datum.c2.x), toString(datum.c2.y), toString(datum.t.x), toString(datum.t.y)], ' ')",
                        type="formula"
                    }
                ]
            },
            {
                name="fLocs",
                source="edges",
                transform=[
                    {
                        expr="datum.f != 0",
                        type="filter"
                    },
                    {
                        fields=[ "s.x", "s.y", "t.x", "t.y", "c1.x", "c1.y", "c2.x", "c2.y" ],
                        ops=[ "mean", "mean", "mean", "mean", "mean", "mean", "mean", "mean" ],
                        as=[ "sx", "sy", "tx", "ty", "c1x", "c1y", "c2x", "c2y" ],
                        groupby=[ "f" ],
                        type="aggregate"
                    }
                ]
            }
        ],
        signals=[
            {
                name="vSize",
                bind={
                    step=10,
                    max=1000,
                    min=100,
                    input="range"
                },
                value=500
            },
            {
                name="strokeWidth",
                bind={
                    step=0.05,
                    max=3,
                    min=0,
                    input="range"
                },
                value=1.5
            },
            {
                name="edgeWidth",
                bind={
                    step=0.05,
                    max=3,
                    min=0,
                    input="range"
                },
                value=1.5
            },
            {
                name="ctRadius",
                bind={
                    step=5,
                    max=100,
                    min=0,
                    input="range"
                },
                value=50
            },
            {
                name="whichSymbol",
                on=[
                    {
                        events="symbol:mousedown",
                        update="datum"
                    },
                    {
                        events="*:mouseup",
                        update="{}"
                    }
                ],
                value={}
            },
            {
                name="newLoc",
                on=[
                    {
                        events="symbol:mouseout[!event.buttons], window:mouseup",
                        update="false"
                    },
                    {
                        events="symbol:mouseover",
                        update="{x: x(), y: y()}"
                    },
                    {
                        events="[symbol:mousedown, window:mouseup] > window:mousemove!",
                        update="{x: x(), y: y()}"
                    }
                ],
                value=false,
                description="State variable for active node newLoc status."
            },
            {
                name="Clicked",
                on=[
                    {
                        events="symbol:mouseout[!event.buttons], window:mouseup",
                        update="false"
                    },
                    {
                        events="symbol:click",
                        update="true"
                    }
                ],
                value=false,
                description="State variable for active node newLoc status."
            }
        ],
        marks=[
            {
                encode={
                    update={
                        strokeWidth={ signal="edgeWidth" },
                        path={ field="mypath" }
                    },
                    enter={
                        stroke={ value="black" }
                    }
                },
                from={ data="edges" },
                type="path"
            },
            {
                encode={
                    update={
                        stroke={ value="black" },
                        x={ field="x" },
                        strokeWidth={ signal="strokeWidth" },
                        size={ signal="vSize" },
                        y={ field="y" }
                    },
                    enter={
                        shape={ field="shape" },
                        fill={ field="color" }
                    }
                },
                from={ data="vertices" },
                type="symbol"
            },
            {
                encode={
                    update={
                        x={ signal="datum.sx*1/2 + datum.tx*1/2" },
                        text={ field="f" },
                        y={ signal="datum.sy*1/2 + datum.ty*1/2" }
                    },
                    enter={
                        baseline={ value="middle" },
                        align={ value="center" }
                    }
                },
                from={ data="fLocs" },
                type="text"
            },
            {
                encode={
                    update={
                        x={ field="m.x" },
                        y={ field="m.y" }
                    }
                },
                from={ data="edges" },
                type="symbol"
            }
        ]
    )
end
