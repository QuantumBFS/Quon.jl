struct Match{R, P}
    parent::Tait{P}
    vertices::Vector{Int}
    half_edges::Vector{Int} # half edge id
end

function Match{R}(parent::Tait{P}, vertices, half_edges) where {R, P}
    Match{R, P}(parent, vertices, half_edges)
end

function Base.show(io::IO, m::Match)
    print(io, "Match:")
    if !isempty(m.vertices)
        println(io)
        print(io, "  vertices = ", m.vertices)
    end

    if !isempty(m.half_edges)
        println(io)
        print(io, "  half_edge = ", m.half_edges)
    end
    return
end

struct Rule{T} end
Rule(r::Symbol) = Rule{r}()

Base.match(r::Rule{R}, tait::Tait{P}) where {R, P} = match!(Match{R, P}[], r, tait)

function match!(matches, ::Rule{:string_genus}, tait::Tait)
    for (v, _) in tait.g.vs_isolated
        v in tait.genuses && push!(matches, Match{:string_genus}(tait, [v], Int[]))
    end
    return matches
end

function match!(matches, ::Rule{:yang_baxter_star}, tait::Tait)
    for (v, _) in tait.g.v2he
        is_open(tait, v) && continue
        hes = trace_vertex(tait, v)
        if length(hes) == 3 && all(x->!is_open(tait, dst(tait, x)), hes)
            push!(matches, Match{:yang_baxter_star}(tait, [v], hes))
        end
    end
    return matches
end

function match!(matches, ::Rule{:yang_baxter_triangle}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        length(hes) == 3 || continue
        has_open_half_edge(tait, hes) || continue
        push!(matches, Match{:yang_baxter_triangle}(tait, [], hes))
    end
    return matches
end

function match!(matches, ::Rule{:charge_rm_v}, tait::Tait)
    for (v, _) in tait.g.v2he
        is_open(tait, v) && continue
        hes = trace_vertex(tait, v)
        ismatch = true
        for he in hes
            haskey(tait.phases, he) || (ismatch = false; break)
            p = phase(tait, he)
            if !p.isparallel
                if !is_singular_change_direction(p)
                    change_direction!(tait, he)
                else
                    ismatch = false; break
                end
            end
            isapprox(real(p.param), 0; atol = quon_atol) || (ismatch = false; break)
            isapprox(rem2pi(imag(p.param), RoundDown), π; atol = quon_atol) || (ismatch = false; break)
        end
        ismatch && push!(matches, Match{:charge_rm_v}(tait, [v], hes))
    end
    return matches
end

function match!(matches, ::Rule{:charge_rm_f}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        has_open_half_edge(tait, hes) || continue
        ismatch = true
        for he in hes
            haskey(tait.phases, he) || (ismatch = false; break)
            p = phase(tait, he)
            if p.isparallel
                if !is_singular_change_direction(p)
                    change_direction!(tait, he)
                else
                    ismatch = false; break
                end
            end
            isapprox(real(p.param), 0; atol = quon_atol) || (ismatch = false; break)
            isapprox(rem2pi(imag(p.param), RoundDown), π; atol = quon_atol) || (ismatch = false; break)
        end
        push!(matches, Match{:charge_rm_f}(tait, [], hes))
    end
    return matches
end

function match!(matches, ::Rule{:z_fusion}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        length(hes) == 2 || continue
        has_open_half_edge(tait, hes) || continue
        push!(matches, Match{:z_fusion}(tait, [], hes))
    end
    return matches
end

function match!(matches, ::Rule{:x_fusion}, tait::Tait)
    for v in vertices(tait)
        is_genus(tait, v) && continue
        hes = trace_vertex(tait, v)
        length(hes) == 2 || continue
        has_open_half_edge(tait, hes) || continue
        push!(matches, Match{:x_fusion}(tait, [], hes))
    end
    return matches
end

# TODO: can we just match 3 genuses?
function match!(matches, ::Rule{:perm_rz}, tait::Tait)
    for v in vertices(tait)
        is_open(tait, v) && continue
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
    return matches
end

function match!(matches, ::Rule{:identity}, tait::Tait)
    for (he_id, theta) in tait.phases
        if is_phase_approx_zero(theta)
            push!(matches, Match{:identity}(tait, [], [he_id]))
        end
    end
    return matches
end

function match!(matches, ::Rule{:genus_fusion}, tait::Tait)
    for g in genuses(tait)
        is_isolated(tait, g) && continue
        hes_g = trace_vertex(tait, g)
        for he_g in hes_g
            face(tait, he_g) == 0 && continue
            he = next(tait, he_g)
            while he != he_g
                v = src(tait, he)
                if is_genus(tait, v)
                    push!(matches, Match{:genus_fusion}(tait, [g, v], [he_g, he]))
                end
                he = next(tait, he)
            end
        end
    end
    return matches
end

function has_open_half_edge(tait::Tait, hes)
    all(hes) do he
        !is_open(tait, dst(tait, he)) && !is_open(tait, src(tait, he))
    end
end

function is_phase_approx_zero(theta::Phase)
    return isapprox(0, theta.param; atol = quon_atol) || isapprox(0, change_direction(theta.param); atol = quon_atol)
end
