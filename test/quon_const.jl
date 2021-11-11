using Test
using Quon: QuonConst, Zero, Pi, HalfPi, NegHalfPi, InfZero, InfPi
using Quon
using Yao

to_Rz(p::Complex) = Rz(p/im)
function to_Rz(p::QuonConst)
    p === Pi && return Rz(π)
    p === Zero && return Rz(0.0)
    p === HalfPi && return Rz(π/2)
    p === NegHalfPi && return Rz(-π/2)
    p === InfZero && return matblock(projector(0))
    p === InfPi && return matblock(projector(1))
end
to_Rx(p::Complex) = Rx(p/im)
function to_Rx(p::QuonConst)
    p === Pi && return Rx(π)
    p === Zero && return Rx(0.0)
    p === HalfPi && return Rx(π/2)
    p === NegHalfPi && return Rx(-π/2)
    p === InfZero && return matblock((1/sqrt(2) + 0.0im)*[1 1; 1 1])
    p === InfPi && return matblock((1/sqrt(2) + 0.0im)*[1 -1; -1 1])
end

function test_quon_const(a, b, c)
    circ1 = chain(1, to_Rz(a), to_Rx(b), to_Rz(c))
    p, q, r = yang_baxter_param(a, b, c)
    circ2 = chain(1, to_Rx(p), to_Rz(q), to_Rx(r))
    return operator_fidelity(matblock(mat(circ1)), matblock(mat(circ2)))/
        sqrt(operator_fidelity(matblock(mat(circ1)), matblock(mat(circ1))) * operator_fidelity(matblock(mat(circ2)), matblock(mat(circ2))))
end

AllConsts = (Zero, Pi, HalfPi, NegHalfPi, InfZero, InfPi)
FiniteConsts = (Zero, Pi, HalfPi, NegHalfPi)

test_quon_const(Pi, Pi, Pi)

test_123 = [test_quon_const(p1, p2, p3) for p1 in AllConsts, p2 in AllConsts, p3 in AllConsts]
test_12 = [test_quon_const(p1, p2, Complex(p3)) for p1 in AllConsts, p2 in AllConsts, p3 in FiniteConsts]
test_23 = [test_quon_const(Complex(p1), p2, p3) for p1 in FiniteConsts, p2 in AllConsts, p3 in AllConsts]
test_13 = [test_quon_const(p1, Complex(p2), p3) for p1 in AllConsts, p2 in FiniteConsts, p3 in AllConsts]
test_1 = [test_quon_const(p1, Complex(p2), Complex(p3)) for p1 in AllConsts, p2 in FiniteConsts, p3 in FiniteConsts]
test_2 = [test_quon_const(Complex(p1), p2, Complex(p3)) for p1 in FiniteConsts, p2 in AllConsts, p3 in FiniteConsts]
test_3 = [test_quon_const(Complex(p1), Complex(p2), p3) for p1 in FiniteConsts, p2 in FiniteConsts, p3 in AllConsts]

@test sum(test_123 .≈ 1) == length(test_123) - 4
@test sum(test_12 .≈ 1) == length(test_12)
@test sum(test_13 .≈ 1) == length(test_13) - 4
@test sum(test_23 .≈ 1) == length(test_23)
@test sum(test_1 .≈ 1) == length(test_1)
@test sum(test_2 .≈ 1) == length(test_2)
@test sum(test_3 .≈ 1) == length(test_3)
