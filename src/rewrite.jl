const rewrite_rules = [
    :string_genus,
    :yang_baxter_star,
    :yang_baxter_triangle,
    :charge_rm_v,
    :charge_rm_f,
    :z_fusion,
    :x_fusion,
    :perm_rz,
    :identity,
    :genus_fusion,
    :swap_genus,
]

function rewrite!(tait::Tait, m::Match{:string_genus})
    rem_vertex!(tait, m.vertices[1])
    return tait
end

function rewrite!(tait::Tait, m::Match{:yang_baxter_star})
    hes = m.half_edges
    v0 = m.vertices[1]
    vs = [dst(tait, he) for he in hes]
    p1, p2, p3 = (quon_param(tait, he) for he in hes)
    q1, q2, q3 = yang_baxter_param_inv(p1, p2, p3)
    add_edge!(tait, vs[1], vs[2], face(tait, hes[1]), q1)
    add_edge!(tait, vs[1], vs[3], face(tait, hes[3]), q2)
    add_edge!(tait, vs[2], vs[3], face(tait, hes[2]), q3)
    rem_vertex!(tait, v0; update = true)
    return tait
end

function rewrite!(tait::Tait, m::Match{:yang_baxter_triangle})
    hes = m.half_edges
    vs = [src(tait, he) for he in hes]
    f = face(tait, hes[1])
    for ids in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
        p1, p2, p3 = (quon_param(tait, hes[ids[1]]), quon_param(tait, hes[ids[3]]), quon_param(tait, hes[ids[2]]))
        if !is_singular_yang_baxter(p1, p2, p3)
            q1, q2, q3 = yang_baxter_param(p1, p2, p3)
            new_v = add_vertex!(tait, f)
            add_edge!(tait, vs[ids[1]], new_v, face(tait, hes[ids[1]]), q1)
            add_edge!(tait, vs[ids[2]], new_v, face(tait, hes[ids[2]]), q2)
            add_edge!(tait, vs[ids[3]], new_v, face(tait, hes[ids[3]]), q3)
            for he in hes 
                rem_edge!(tait, he; update = true)
            end
            break
        end
    end
    return tait
end

function rewrite!(tait::Tait{QuonParam{T}}, m::Match{:charge_rm_v}) where T
    hes = m.half_edges
    v = m.vertices[1]
    p0 = QuonParam{T}(zero(T), true)
    if (σ_inv(tait, hes[1]) == hes[end]) || (length(hes) == 1)
        for he in hes
            tait.quon_params[he] = p0
            tait.quon_params[twin(tait, he)] = p0
            rewrite!(tait, Match{:identity}([], [he]))
        end
        return tait
    end

    hes_identity = Int[]
    for i = 1:(length(hes)-1)
        he1 = hes[i]
        he2 = hes[i+1]
        v1 = dst(tait, he1)
        v2 = dst(tait, he2)
        f = face(tait, he1)
        push!(hes_identity, add_edge!(tait, v1, v2, f, p0)[1])
    end
    for i = 1:(length(hes))
        rem_edge!(tait, hes[i]; update = true)
        rewrite!(tait, Match{:identity}([], [hes_identity[i]]))
    end
    return tait
end

function rewrite!(tait::Tait{QuonParam{T}}, m::Match{:charge_rm_f}) where T
    hes = m.half_edges
    f = face(tait, hes[1])
    if next(tait, hes[end]) != hes[1]
        add_edge!(tait, src(tait, hes[1]), dst(tait, hes[end]), f, QuonParam{T}(Pi, false))
    end
    for he in hes
        rem_edge!(tait, he; update = true)
    end
    return tait
end

function rewrite!(tait::Tait{QuonParam{T}}, m::Match{:z_fusion}) where T
    he1, he2 = m.half_edges
    twin1 = twin(tait, he1)
    twin2 = twin(tait, he2)
    p1 = quon_param(tait, he1)
    p2 = quon_param(tait, he2)
    p = add_orthorgonal(p1, p2)
    p0 = QuonParam{T}(zero(p2.param), false)
    tait.quon_params[he1] = p
    tait.quon_params[twin1] = p
    tait.quon_params[he2] = p0
    tait.quon_params[twin2] = p0
    rewrite!(tait, Match{:identity}([], [he2]))
    return tait
end

function rewrite!(tait::Tait{QuonParam{T}}, m::Match{:x_fusion}) where T
    he1, he2 = m.half_edges
    twin1 = twin(tait, he1)
    twin2 = twin(tait, he2)
    p1 = quon_param(tait, he1)
    p2 = quon_param(tait, he2)
    p = add_parallel(p1, p2)
    p0 = QuonParam{T}(zero(p2.param), true)
    tait.quon_params[he1] = p
    tait.quon_params[twin1] = p
    tait.quon_params[he2] = p0
    tait.quon_params[twin2] = p0
    rewrite!(tait, Match{:identity}([], [he2]))
    return tait
end

function rewrite!(tait::Tait, m::Match{:perm_rz})
    v_ct, v_g = m.vertices
    f = face(tait, m.half_edges[2])
    old_he = m.half_edges[1]
    add_edge!(tait, v_ct, v_g, f, quon_param(tait, old_he))
    rem_edge!(tait, old_he; update = true)
    return tait
end

function rewrite!(tait::Tait, m::Match{:identity})
    he = m.half_edges[1]
    p = quon_param(tait, he)
    is_zero(p.param) || (p = change_direction(p))
    if is_parallel(p)
        contract_edge!(tait, he)
    else
        rem_edge!(tait, he; update = true)
    end
    return tait
end

function rewrite!(tait::Tait{QuonParam{T}}, m::Match{:genus_fusion}) where T
    he_g1, _ = m.half_edges
    g1, g2 = m.vertices
    f = face(tait, he_g1)
    new_he1, _ = add_edge!(tait, g1, g2, f, QuonParam{T}(zero(T), true))
    rewrite!(tait, Match{:identity}([], [new_he1]))
end

function rewrite!(tait::Tait, m::Match{:swap_genus})
    g = m.vertices[1]
    delete!(tait.genuses, g)
    return tait
end