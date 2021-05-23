function rewrite!(tait::Tait, m::Match{:string_genus})
    rem_vertex!(tait, m.vertices[1])
    return tait
end

function rewrite!(tait::Tait, m::Match{:yang_baxter_star})
    hes = m.half_edges
    v0 = m.vertices[1]
    vs = [dst(tait, he) for he in hes]
    fs = [face(tait, he) for he in hes]
    p1, p2, p3 = (phase(tait, he) for he in hes)
    q1, q2, q3 = yang_baxter_param_inv(p1, p2, p3)
    add_edge!(tait, vs[1], vs[2], fs[1], q1)
    add_edge!(tait, vs[1], vs[3], fs[3], q2)
    add_edge!(tait, vs[2], vs[3], fs[2], q3)
    rem_vertex!(tait, v0; update = true)
    return tait
end

function rewrite!(tait::Tait, m::Match{:yang_baxter_triangle})
    hes = m.half_edges
    vs = [src(tait, he) for he in hes]
    f = face(tait, hes[1])
    new_v = add_vertex!(tait, f)
    p1, p2, p3 = (phase(tait, hes[1]), phase(tait, hes[3]), phase(tait, hes[2]))
    q1, q2, q3 = yang_baxter_param(p1, p2, p3)
    add_edge!(tait, vs[1], new_v, face(tait, hes[1]), q1)
    add_edge!(tait, vs[2], new_v, face(tait, hes[2]), q2)
    add_edge!(tait, vs[3], new_v, face(tait, hes[3]), q3)
    for he in hes 
        rem_edge!(tait, he; update = true)
    end
    return tait
end

function rewrite!(tait::Tait, m::Match{:z_fusion})
    he1, he2 = m.half_edges
    twin1 = twin(tait, he1)
    twin2 = twin(tait, he2)
    p1 = phase(tait, he1)
    p2 = phase(tait, he2)
    p = p1 + p2
    p0 = Phase{typeof(p2.param)}(zero(p2.param), p2.isparallel)
    tait.phases[he1] = p
    tait.phases[twin1] = p
    tait.phases[he2] = p0
    tait.phases[twin2] = p0
    rewrite!(tait, Match{:identity}(tait, [], [he2]))
    return tait
end

rewrite!(tait::Tait, m::Match{:x_fusion}) = rewrite!(tait, Match{:z_fusion}(m.parent, m.vertices, m.half_edges))

function rewrite!(tait::Tait, m::Match{:perm_rz})
    v_ct, v_g = m.vertices
    f = face(tait, m.half_edges[2])
    old_he = m.half_edges[1]
    add_edge!(tait, v_ct, v_g, f, phase(tait, old_he))
    rem_edge!(tait, old_he)
    return tait
end

function rewrite!(tait::Tait, m::Match{:identity})
    he = m.half_edges[1]
    p = phase(tait, he)
    iszero(p.param) || (p = change_direction(p))
    if p.isparallel
        contract_edge!(tait, he)
    else
        rem_edge!(tait, he; update = true)
    end
    return tait
end
