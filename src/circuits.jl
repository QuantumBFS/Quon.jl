function Tait(n::Integer)
    circ = tait_id()
    id = tait_id()
    for _ in 1:(n-1)
        tensor_product!(circ, id)
    end
    return circ
end

using ZXCalculus: push_gate!

function ZXCalculus.push_gate!(q::Tait, ::Val{:H}, loc)
    h = tait_hadamard()
    contract!(q, h, [q.outputs[loc]], copy(h.inputs))
    return q
end

function ZXCalculus.push_gate!(q::Tait, ::Val{:Rz}, loc, p)
    g = tait_rz(p*im)
    contract!(q, g, [q.outputs[loc]], copy(g.inputs))
    return q
end
ZXCalculus.push_gate!(q::Tait, ::Val{:shift}, loc, p) = ZXCalculus.push_gate!(q, Val(:Rz), loc, p)

function ZXCalculus.push_gate!(q::Tait, ::Val{:Rx}, loc, p)
    g = tait_rx(p*im)
    contract!(q, g, [q.outputs[loc]], copy(g.inputs))
    return q
end

function ZXCalculus.push_gate!(q::Tait, ::Val{:SWAP}, loc1, loc)
    g = tait_swap()
    # contract!(q, g, [q.])
end

function ZXCalculus.push_gate!(q::Tait, ::Val{:CNOT}, loc, ctrl)
    g = tait_cnot(ctrl, loc)
    contract!(q, g, q.outputs[min(loc, ctrl):max(loc, ctrl)], copy(g.inputs))
    return q
end

function ZXCalculus.push_gate!(q::Tait, ::Val{:CZ}, loc, ctrl)
    g = tait_cz(ctrl, loc)
    contract!(q, g, q.outputs[min(loc, ctrl):max(loc, ctrl)], copy(g.inputs))
    return q
end
