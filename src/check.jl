function has_subgraph(tait::Tait, m::Match)
    all(has_vertex(tait, v) for v in m.vertices) || return false
    all(has_half_edge(tait, he) for he in m.half_edges) || return false
    return true
end

check(tait::Tait, m::Match) = false

function check(tait::Tait, m::Match{:string_genus})
    v = m.vertices[1]
    return v in tait.genuses
end

function check(tait::Tait, m::Match{:yang_baxter_star})
    has_subgraph(tait, m) || return false
    v = m.vertices[1]

    (is_open_vertex(tait, v) || is_genus(tait, v)) && return false
    hes = trace_vertex(tait, v)
    if length(hes) == 3 
        nbs = Set(dst(tait, he) for he in hes)
        if length(nbs) == 3 && all(x->!is_open_vertex(tait, dst(tait, x)), hes)
            return Set(hes) == Set(m.half_edges)
        end
    end
    return false
end

function check(tait::Tait, m::Match{:yang_baxter_triangle})
    has_subgraph(tait, m) || return false
    hes = m.half_edges

    length(hes) == 3 || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false

    f = face(tait, hes[1])
    f != 0 || return false
    f == face(tait, hes[2]) == face(tait, hes[3]) || return false
    return length(trace_face(tait, f)) == 3
end

function check(tait::Tait, m::Match{:charge_rm_v}) 
    has_subgraph(tait, m) || return false
    hes = m.half_edges
    v = m.vertices[1]

    is_open_vertex(tait, v) && return false
    src(tait, hes[1]) == v || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false
    for i = lastindex(hes):2
        σ_inv(tait, hes[i]) == hes[i-1] || return false
        p = quon_param(tait, hes[i])
        is_para_pi(p) || return false
    end
    p1 = quon_param(tait, hes[1])
    return is_para_pi(p1)
end

function check(tait::Tait, m::Match{:charge_rm_f})
    has_subgraph(tait, m) || return false
    hes = m.half_edges
    
    face(tait, hes[1]) != 0 || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false
    for i = 1:(lastindex(hes)-1)
        next(tait, hes[i]) == hes[i+1] || return false
        p = quon_param(tait, hes[i])
        is_orth_pi(p) || return false
    end
    p_end = quon_param(tait, hes[end])
    return is_orth_pi(p_end)
end

function check(tait::Tait, m::Match{:z_fusion})
    has_subgraph(tait, m) || return false
    hes = m.half_edges

    length(hes) == 2 || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false
    next(tait, hes[1]) == hes[2] || return false
    next(tait, hes[2]) == hes[1] || return false
    return face(tait, hes[1]) != 0
end

function check(tait::Tait, m::Match{:x_fusion})
    has_subgraph(tait, m) || return false
    v = m.vertices[1]
    hes = m.half_edges
    
    is_genus(tait, v) && return false
    length(hes) == 2 || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false
    src(tait, hes[1]) == v || return false
    σ_inv(tait, hes[1]) == hes[2] || return false
    return σ_inv(tait, hes[2]) == hes[1]
end

function check(tait::Tait, m::Match{:perm_rz})
    has_subgraph(tait, m) || return false
    (v, g2) = m.vertices
    he1, he2 = m.half_edges

    g1 = dst(tait, he1)
    is_open_half_edge(tait, he1) && return false
    src(tait, he2) == g2 || return false
    (is_genus(tait, g1) && is_genus(tait, g2)) || return false
    is_genus(tait, v) && return false
    hes_v = trace_vertex(tait, v)
    fs_v = Set(face(tait, he) for he in hes_v)
    return face(tait, he2) in fs_v
end

function check(tait::Tait, m::Match{:identity})
    has_subgraph(tait, m) || return false
    he = m.half_edges[1]
    return is_para_zero(quon_param(tait, he)) || is_orth_zero(quon_param(tait, he))
end

function check(tait::Tait, m::Match{:genus_fusion})
    has_subgraph(tait, m) || return false
    g1, g2 = m.vertices
    he1, he2 = m.half_edges

    (is_genus(tait, g1) && is_genus(tait, g2)) || return false
    (src(tait, he1) == g1 && src(tait, he2) == g2) || return false
    face(tait, he1) == face(tait, he2) || return false
    return face(tait, he1) != 0
end

function check(tait::Tait, m::Match{:swap_genus})
    has_subgraph(tait, m) || return false
    vs = m.vertices
    hes = m.half_edges

    all(is_genus(tait, g) for g in vs[1:4]) || return false
    all(!is_genus(tait, v) for v in vs[5:9]) || return false
    v0 = vs[9]
    all(src(tait, he) == v0 for he in hes) || return false
    all(length(trace_face(tait, face(tait, he))) == 4 for he in hes) || return false
    all(check_swap_parameters(tait, v) for v in vs[5:9]) || return false
    g1 = vs[1]
    return !is_genus_connected_to_open_edge(tait, g1)
end
