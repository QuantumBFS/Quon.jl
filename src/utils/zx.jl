using Multigraphs: Multigraph
import ZXCalculus
using ZXCalculus: ZXDiagram, SpiderType, add_spider!, Phase

function ZXCalculus.ZXDiagram(nin::Integer, nout::Integer)
    N = nin + nout
    mg = Multigraph(nin + nout)
    st = vcat([SpiderType.In for _ = 1:nin], [SpiderType.Out for _ = 1:nout])
    ps = [Phase(0//1) for _ = 1:N]
    zxd = ZXDiagram(mg, st, ps)
    sort!(zxd.inputs)
    sort!(zxd.outputs)
    return zxd
end

function ZXCalculus.ZXDiagram(q::Tait)
    simplify!(q, Rule(:genus_fusion))
    nin = length(q.inputs)
    nout = length(q.outputs)
    zxd = ZXDiagram(nin, nout)
    v_map = Dict{Int, Int}()
    for i = 1:nin
        v_map[q.inputs[i]] = i
    end
    for i = 1:nout
        v_map[q.outputs[i]] = i + nin
    end

    for v in vertices(q)
        haskey(v_map, v) && continue
        is_genus(q, v) && continue
        v_map[v] = add_spider!(zxd, SpiderType.Z)
    end

    he_map = Dict{Int, Int}()
    for he in half_edges(q)
        he_twin = twin(q, he)
        s = src(q, he)
        d = dst(q, he)
        if !haskey(he_map, he)
            if haskey(q.quon_params, he)
                p = quon_param(q, he)
                if is_genus(q, s) || is_genus(q, d)
                    v_z = v_map[is_genus(q, s) ? d : s]
                    if real(p.param) == 0
                        if !is_parallel(p)
                            zxd.ps[v_z] += imag(p.param) / pi
                            he_map[he] = v_z
                            he_map[he_twin] = v_z
                        else
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p.param) / pi), [v_z])
                            he_map[he] = v_x
                            he_map[he_twin] = v_x
                        end
                    else    # the parameter is not pure imag number
                        p1 = QuonParam{QuonComplex}(real(p.param) + pi/2*im, is_parallel(p))
                        p2 = QuonParam{QuonComplex}(imag(p.param) - pi/2*im, is_parallel(p))
                        p1 = change_direction(p1)
                        if !is_parallel(p1)
                            zxd.ps[v_z] += imag(p1.param) / pi
                            he_map[he] = v_z
                            he_map[he_twin] = v_z
                        else
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p1.param) / pi), [v_z])
                            he_map[he] = v_x
                            he_map[he_twin] = v_x
                        end
                        if !is_parallel(p2)
                            zxd.ps[v_z] += imag(p2.param) / pi
                            he_map[he] = v_z
                            he_map[he_twin] = v_z
                        else
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p2.param) / pi), [v_z])
                            he_map[he] = v_x
                            he_map[he_twin] = v_x
                        end
                    end
                else
                    if real(p.param) == 0
                        if is_parallel(p)
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p.param) / pi), [v_map[s], v_map[d]])
                            he_map[he] = v_x
                            he_map[he_twin] = v_x
                        else
                            v_x = add_spider!(zxd, SpiderType.X, Phase(0//1), [v_map[s], v_map[d]])
                            v_z = add_spider!(zxd, SpiderType.Z, Phase(imag(p.param) / pi), [v_x])
                            he_map[he] = v_z
                            he_map[he_twin] = v_z
                        end
                    else    # the parameter is not pure imag number
                        p1 = QuonParam{QuonComplex}(real(p.param) + pi/2*im, is_parallel(p))
                        p2 = QuonParam{QuonComplex}(imag(p.param) - pi/2*im, is_parallel(p))
                        p1 = change_direction(p1)
                        if is_parallel(p1)
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p1.param) / pi), [v_map[s], v_map[d]])
                            v_z = add_spider!(zxd, SpiderType.Z, Phase(imag(p2.param) / pi), [v_x])
                            he_map[he] = v_x
                            he_map[he_twin] = v_x
                        else
                            v_x = add_spider!(zxd, SpiderType.X, Phase(imag(p2.param) / pi), [v_map[s], v_map[d]])
                            v_z = add_spider!(zxd, SpiderType.Z, Phase(imag(p1.param) / pi), [v_x])
                            he_map[he] = v_z
                            he_map[he_twin] = v_z
                        end
                    end
                end
            else
                if !(is_genus(q, s) || is_genus(q, d))
                    ZXCalculus.add_edge!(zxd, v_map[s], v_map[d])
                    he_map[he] = v_map[s]
                    he_map[he_twin] = v_map[d]
                end
            end
        end
    end
    return zxd
end