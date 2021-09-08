struct Match{R}
    # TODO: do we need the `parent`
    # parent::Tait{P}
    vertices::Vector{Int}
    half_edges::Vector{Int} # half edge id
end

# function Match{R}(parent::Tait{P}, vertices, half_edges) where {R, P}
#     Match{R, P}(parent, vertices, half_edges)
# end

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

"""
    Rule{R}

A struct for rules. `R` can be `:string_genus`, `yang_baxter_star`, `yang_baxter_triangle`, 
`charge_rm_v`, `charge_rm_f`, `z_fusion`, `x_fusion`, `perm_rz`, `identity`, `genus_fusion`, 
`swap_genus`.
"""
struct Rule{T} end
Rule(r::Symbol) = Rule{r}()

Base.match(r::Rule{R}, tait::Tait{P}) where {R, P} = match!(Match{R}[], r, tait)

function match!(matches, ::Rule{:string_genus}, tait::Tait)
    for (v, _) in tait.g.vs_isolated
        v in tait.genuses && push!(matches, Match{:string_genus}([v], Int[]))
    end
    return matches
end

function match!(matches, ::Rule{:yang_baxter_star}, tait::Tait)
    for v in vertices(tait)
        is_open_vertex(tait, v) && continue
        hes = trace_vertex(tait, v)
        if length(hes) == 3 && all(x->!is_open_vertex(tait, dst(tait, x)), hes)
            push!(matches, Match{:yang_baxter_star}([v], hes))
        end
    end
    return matches
end

function match!(matches, ::Rule{:yang_baxter_triangle}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        length(hes) == 3 || continue
        all(x->!is_open_vertex(tait, dst(tait, x)), hes) || continue
        push!(matches, Match{:yang_baxter_triangle}([], hes))
    end
    return matches
end

function match!(matches, ::Rule{:charge_rm_v}, tait::Tait)
    for v in vertices(tait)
        is_open_vertex(tait, v) && continue
        hes = trace_vertex(tait, v)
        if length(hes) == 1
            if !is_open_half_edge(tait, hes[1])
                p = phase(tait, hes[1])
                if is_phase_pi(p) && is_parallel(p)
                    push!(matches, Match{:charge_rm_v}([v], hes))
                end
            end
            continue
        end
        
        i1 = findfirst(hes) do x
            if is_open_half_edge(tait, x)
                return true
            else 
                return !is_phase_pi(phase(tait, x)) || !is_parallel(phase(tait, x))
            end
        end
        if i1 === nothing
            push!(matches, Match{:charge_rm_v}([v], hes))
            continue
        end

        hes = permute!(hes, [collect(i1:length(hes)); collect(1:(i1-1))])
        hes_match = Int[]
        for i = eachindex(hes)
            is_match = false
            he = hes[i]
            if !is_open_half_edge(tait, he)
                p = phase(tait, he)
                if is_parallel(p) && is_phase_pi(p)
                    push!(hes_match, he)
                    is_match = true
                    if he == hes[end]
                        (length(hes_match) > 1) && push!(matches, Match{:charge_rm_v}([v], hes_match))
                    end
                end
            end
            if !is_match
                if length(hes_match) > 1
                    push!(matches, Match{:charge_rm_v}([v], hes_match))
                end
                !isempty(hes_match) && (hes_match = Int[])
            end
        end
    end
    return matches
end

function match!(matches, ::Rule{:charge_rm_f}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        
        i1 = findfirst(hes) do x
            if is_open_half_edge(tait, x)
                return true
            else 
                return !is_phase_pi(phase(tait, x)) || !is_parallel(phase(tait, x))
            end
        end
        if i1 === nothing
            println(f)
            (length(hes) > 1) && push!(matches, Match{:charge_rm_f}([], hes))
            continue
        end

        hes = permute!(hes, [collect(i1:length(hes)); collect(1:(i1-1))])
        hes_match = Int[]
        for i = eachindex(hes)
            is_match = false
            he = hes[i]
            if !is_open_half_edge(tait, he)
                p = phase(tait, he)
                if !is_parallel(p) && is_phase_pi(p)
                    push!(hes_match, he)
                    is_match = true
                    if he == hes[end]
                        (length(hes_match) > 1) && push!(matches, Match{:charge_rm_f}([], hes_match))
                    end
                end
            end
            if !is_match
                if length(hes_match) > 1
                    push!(matches, Match{:charge_rm_f}([], hes_match))
                end
                !isempty(hes_match) && (hes_match = Int[])
            end
        end
    end
    return matches
end

function match!(matches, ::Rule{:z_fusion}, tait::Tait)
    for f in faces(tait)
        f == 0 && continue
        hes = trace_face(tait, f)
        length(hes) == 2 || continue
        all(he -> !is_open_half_edge(tait, he), hes) || continue
        push!(matches, Match{:z_fusion}([], hes))
    end
    return matches
end

function match!(matches, ::Rule{:x_fusion}, tait::Tait)
    for v in vertices(tait)
        is_genus(tait, v) && continue
        is_open_vertex(tait, v) && continue
        hes = trace_vertex(tait, v)
        length(hes) == 2 || continue
        all(he -> !is_open_half_edge(tait, he), hes) || continue
        push!(matches, Match{:x_fusion}([v], hes))
    end
    return matches
end

function match!(matches, ::Rule{:perm_rz}, tait::Tait)
    for v in vertices(tait)
        is_open_vertex(tait, v) && continue
        hes = trace_vertex(tait, v)
        # length(hes) >= 3 || continue
        ids = findall(hes) do he
            is_genus(tait, dst(tait, he))
        end
        # TODO: Should we check isolated genuses?
        length(ids) == 0 && continue
        for idx in ids
            he_genus = hes[idx]
            # find half_edge on the other face that has a genus
            he_out = σ_inv(tait, he_genus)
            for _ in 2:(length(hes) - 1)
                he_out = σ_inv(tait, he_out)
                he = he_out
                while !is_genus(tait, src(tait, he))
                    he = next(tait, he)
                    he == he_out && break
                end
                if he != he_out && src(tait, he) != dst(tait, he_genus)
                    push!(matches, Match{:perm_rz}([v, src(tait, he)], [he_genus, he]))
                end
            end
        end
    end
    return matches
end

function match!(matches, ::Rule{:identity}, tait::Tait)
    for (he_id, theta) in tait.phases
        if is_phase_zero(theta)
            push!(matches, Match{:identity}([], [he_id]))
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
                    push!(matches, Match{:genus_fusion}([g, v], [he_g, he]))
                end
                he = next(tait, he)
            end
        end
    end
    return matches
end

function match!(matches, ::Rule{:swap_genus}, tait::Tait)
    for v0 in vertices(tait)
        is_genus(tait, v0) && continue
        hes_v0 = trace_vertex(tait, v0)
        length(hes_v0) == 4 || continue
        vs = [dst(tait, he) for he in hes_v0]

        faces_contain_4_edges = true
        for he in hes_v0
            f = face(tait, he)
            f_hes = trace_face(tait, f)
            if length(f_hes) != 4
                faces_contain_4_edges = false
                break
            end
        end
        faces_contain_4_edges || continue

        vertices_contain_4_neighbors = true
        for v in vs
            v_hes = trace_vertex(tait, v)
            if is_genus(tait, v) || length(v_hes) != 4
                vertices_contain_4_neighbors = false
                break
            end
        end
        vertices_contain_4_neighbors || continue

        # check parameters
        check_swap_parameters(tait, v0) || continue
        for he in hes_v0
            v = dst(tait, he)
            check_swap_parameters(tait, v) || continue
        end

        has_4_connected_genuses = true
        connected_genuses = Int[]
        for he in hes_v0
            he_twin = twin(tait, he)
            he_g = σ_inv(tait, σ_inv(tait, he_twin))
            v_g = dst(tait, he_g)
            if is_genus(tait, v_g)
                push!(connected_genuses, v_g)
            else
                has_4_connected_genuses = false
                break
            end
        end
        has_4_connected_genuses || continue

        # find genus has open edge
        idx = findfirst(x->!is_genus_connected_to_open_edge(tait, x), connected_genuses)
        idx === nothing && continue

        push!(matches, 
            Match{:swap_genus}(
                [connected_genuses[idx:end]; connected_genuses[1:idx-1]; 
                    vs[idx:end]; vs[1:idx-1]; v0
                ], 
                [hes_v0[idx:end]; hes_v0[1:idx-1]]
            )
        )        
    end
    return matches
end

function check_swap_parameters(tait::Tait, v)
    hes = trace_vertex(tait, v)
    length(hes) == 4 || return false
    any(is_open_half_edge(tait, he) for he in hes) && return false
    (p1, p2, p3, p4) = [phase(tait, he) for he in hes]
    if is_phase_para_half_pi(p1) && is_phase_para_half_pi(p3) && 
        is_phase_perp_half_pi(p2) && is_phase_perp_half_pi(p4) 
        return true
    elseif is_phase_para_half_pi(p2) && is_phase_para_half_pi(p4) && 
        is_phase_perp_half_pi(p1) && is_phase_perp_half_pi(p3) 
        return true
    end
    return false
end

function is_genus_connected_to_open_edge(tait::Tait, genus)
    for he in trace_vertex(tait, genus)
        if is_open_vertex(tait, dst(tait, he))
            return true
        end
    end
    return false
end
