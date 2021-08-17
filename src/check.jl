check(tait::Tait{P}, m::Match{R, P}) where {R, P} = false

function check(tait::Tait{P}, m::Match{:string_genus, P}) where P
    v = m.vertices[1]
    return v in tait.genuses
end

function check(tait::Tait{P}, m::Match{:yang_baxter_star, P}) where P
    v = m.vertices[1]
    v in vertices(tait) || return false
    is_open(tait, v) && return false
    hes = trace_vertex(tait, v)
    if length(hes) == 3 && all(x->!is_open(tait, dst(tait, x)), hes)
        return true
    end
    return false
end

function check(tait::Tait{P}, m::Match{:yang_baxter_triangle, P}) where P
    hes = m.half_edges
    length(hes) == 3 || return false
    has_open_half_edge(tait, hes) || return false

    face(tait, hes[1]) == face(tait, hes[2]) == face(tait, hes[3]) || return false
    return face(tait, hes[1]) != 0
end

function check(tait::Tait{P}, m::Match{:charge_rm_v, P}) where P 
end

function check(tait::Tait{P}, m::Match{:charge_rm_f, P}) where P
end

function check(tait::Tait{P}, m::Match{:z_fusion, P}) where P
    hes = m.half_edges
    length(hes) == 2 || return false
    has_open_half_edge(tait, hes) || return false
    face(tait, hes[1]) == face(tait, hes[2]) || return false
    return face(tait, hes[1]) != 0
end

function check(tait::Tait{P}, m::Match{:x_fusion, P}) where P
    v = m.vertices[1]
    v in vertices(tait) || return false
    is_genus(tait, v) && return false
    hes = trace_vertex(tait, v)
    length(hes) == 2 || return false
    has_open_half_edge(tait, hes) || return false
    return true
end

function check(tait::Tait{P}, m::Match{:perm_rz, P}) where P
    v = m.vertices[1]
    v in vertices(tait) || return false
    

    is_open(tait, v) && return false
    hes = trace_vertex(tait, v)
    length(hes) >= 3 || continue
    ids = findall(hes) do he
        is_genus(tait, dst(tait, he))
    end
    length(ids) == 0 && continue
    for idx in ids
        he_genus = hes[idx]
        # find half_edge on the other face that has a genus
        he_out = he_genus
        for _ in 2:(length(hes) - 1)
            he_out = σ_inv(tait, he_out)
            he = twin(tait, he_out)
            he0 = he
            while !is_genus(tait, dst(tait, he))
                he = next(tait, he)
                he == he0 && break
            end
            if he != he0 && dst(tait, he) != dst(tait, he_genus)
                push!(matches, Match{:perm_rz}(
                    tait,
                    [v, dst(tait, he)],
                    [he_genus, he]
                    )
                )
            end
        end
    end
end

function check(tait::Tait{P}, m::Match{:identity, P}) where P
end

function check(tait::Tait{P}, m::Match{:genus_fusion, P}) where P
end

function check(tait::Tait{P}, m::Match{:swap_genus, P}) where P
end
