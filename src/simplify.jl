export simplify!

function simulated_annealing(f, x0, temp)
    l_old = f(x0)
    loss = l_old
    @show loss
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
        loss = l_old
        @show loss
        l_old = l_new
    end
    return candidate
end

function simplify_instance(tait) # -> cost
    return function instance(rules::Vector{Symbol}) # -> cost
        for rule in rules
            matches = match(Rule(rule), tait)
            isempty(matches) && continue
            m = first(matches)
            # @show rule, m
            rewrite!(tait, m)
        end
        return tait_loss(tait)
    end
end

function tait_loss(q::Tait)
    return length(q.quon_params)/2 + length(q.genuses) * 10
end

function simplify!(tait::Tait)
    x0 = rand(rewrite_rules, 10)
    instance = simplify_instance(tait)
    rules = simulated_annealing(instance, x0, 1.0)
    return tait
end

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
