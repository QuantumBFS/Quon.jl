struct Match{P}
    parent::Tait{P}
    vertices::Vector{Int}
    half_edges::Vector{Int} # half edge id
end

struct Rule{T} end

match(r::Rule, tait::Tait{P}) where P = match!(Match{P}[], r, tait)

function match!(matches, ::Rule{:string_genus}, tait::Tait)
    for v in tait.genuses

    end
    return matches
end

function match!(matches, ::Rule{:yang_baxter_star}, tait::Tait)
    
    return matches
end

function match!(matches, ::Rule{:yang_baxter_triangle}, tait::Tait)
    return matches
end

function match!(matches, ::Rule{:fusion}, tait::Tait)
    return matches
end

function match!(matches, ::Rule{:identity}, tait::Tait)
    return matches
end
