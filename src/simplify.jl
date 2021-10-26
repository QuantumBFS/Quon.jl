export simplify!


function simulated_annealing(f, x0, temp)
    l_old = f(x0)
    candidate = copy(x0)
    for epoch in 1:1000
        # local update
        idx = rand(1:length(candidate))
        candidate_idx_old = candidate[idx]
        candidate[idx] = rand(rewrite_rules)

        l_new = f(candidate)
        t = temp / epoch
        ratio = exp((l_old - l_new)/t)
        if l_new < l_old || rand() < ratio
            # pick candidate
        else
            candidate[idx] = candidate_idx_old
        end
        l_old = l_new
    end
    return candidate
end

function simplify_instance(tait) # -> cost
    return function instance(rules::Vector{Symbol}) # -> cost
        for rule in rules
            matches = match(Rule(rule), tait)
            isempty(matches) && continue
            rewrite!(tait, first(matches))
        end
        return tait_loss(tait)
    end
end

function tait_loss(q::Tait)
    return length(q.phases)/2 + length(q.genuses) * 10
end

function simplify!(tait::Tait)
    x0 = rand(rewrite_rules, 10)
    instance = simplify_instance(tait)
    rules = simulated_annealing(instance, x0, 1.0)
    return tait
end

# function simplify!(q::Tait, r::Rule{R}) where {R}
#     ms = match(r, q)
#     while length(ms) > 0
#         for m in ms
#             if check(q, m)
#                 rewrite!(q, m)
#             end
#         end
#         ms = match(r, q)
#     end
#     return q
# end

# function simplify!(q::Tait, rs::Vector{Rule})
#     ms = [match(r, q) for r in rs]
#     while any(length.(ms) .> 0)
#         for ms_r in ms
            
#         end
#     end
# end

# function naive_simplify!(q::Tait)
#     simplify!(q, Rule(:z_fusion))
#     simplify!(q, Rule(:x_fusion))
#     simplify!(q, Rule(:identity))

#     m_yb_star = match(Rule(:yang_baxter_star), q)
#     @show m_yb_star
#     if length(m_yb_star) > 0
#         for m in m_yb_star
#             if check(q, m)
#                 can_simp = [length(trace_face(q, face(q, he))) == 3 for he in m.half_edges]
#                 @show can_simp
#                 if any(can_simp)
#                     rewrite!(q, m)
#                 end
#             end
#         end
#     end
#     return q
# end
# export naive_simplify!