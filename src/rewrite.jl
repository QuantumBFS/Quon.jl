function rewrite!(tait::Tait, match::Match{:string_genus})
end

function rewrite!(tait::Tait, match::Match{:yang_baxter_star})
end

function rewrite!(tait::Tait, match::Match{:yang_baxter_triangle})
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
