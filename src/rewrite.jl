function rewrite!(tait::Tait, match::Match{:string_genus})
end

function rewrite!(tait::Tait, match::Match{:yang_baxter_star})
    hes = match.half_edges
    v0 = match.vertices[1]
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

function rewrite!(tait::Tait, match::Match{:yang_baxter_triangle})
    hes = match.half_edges
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

function rewrite!(tait::Tait, match::Match{:z_fusion})
end

function rewrite!(tait::Tait, match::Match{:x_fusion})
end

function rewrite!(tait::Tait, match::Match{:perm_rz})
    v_ct, v_g = match.vertices
    f = face(tait, match.half_edges[2])
    old_he = match.half_edges[1]
    add_edge!(tait, v_ct, v_g, f, phase(tait, old_he))
    rem_edge!(tait, old_he)
    return tait
end

function rewrite!(tait::Tait, match::Match{:identity})
end
