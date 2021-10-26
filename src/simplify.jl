function simplify!(q::Tait, r::Rule{R}) where {R}
    ms = match(r, q)
    while length(ms) > 0
        for m in ms
            if check(q, m)
                rewrite!(q, m)
            end
        end
        ms = match(r, q)
    end
    return q
end

function simplify!(q::Tait, rs::Vector{Rule})
    ms = [match(r, q) for r in rs]
    while any(length.(ms) .> 0)
        for ms_r in ms
            
        end
    end
end

function naive_simplify!(q::Tait)
    simplify!(q, Rule(:z_fusion))
    simplify!(q, Rule(:x_fusion))
    simplify!(q, Rule(:identity))

    m_yb_star = match(Rule(:yang_baxter_star), q)
    @show m_yb_star
    if length(m_yb_star) > 0
        for m in m_yb_star
            if check(q, m)
                can_simp = [length(trace_face(q, face(q, he))) == 3 for he in m.half_edges]
                @show can_simp
                if any(can_simp)
                    rewrite!(q, m)
                end
            end
        end
    end
    return q
end
export naive_simplify!