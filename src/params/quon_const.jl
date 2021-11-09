struct QuonConst{S} <: AbstractQuonParam end

const ParaHalfPi = QuonConst{Symbol("∥ π/2⋅im")}()
const PerpHalfPi = QuonConst{Symbol("⊥ π/2⋅im")}()
const ParaNegHalfPi = QuonConst{Symbol("∥ -π/2⋅im")}()
const PerpNegHalfPi = QuonConst{Symbol("⊥ -π/2⋅im")}()
const ParaZero = QuonConst{Symbol("∥ 0")}()
const PerpZero = QuonConst{Symbol("⊥ 0")}()
const ParaPi = QuonConst{Symbol("∥ π⋅im")}()
const PerpPi = QuonConst{Symbol("⊥ π⋅im")}()
const ParaInfZero = QuonConst{Symbol("CD(⊥ 0)")}()
const PerpInfZero = QuonConst{Symbol("CD(∥ 0)")}()
const ParaInfPi = QuonConst{Symbol("CD(⊥ π⋅im)")}()
const PerpInfPi = QuonConst{Symbol("CD(∥ π⋅im)")}()

Base.show(io::IO, ::QuonConst{S}) where S = print(io, "$S")
Base.copy(p::QuonConst{S}) where S = p
Base.:(==)(p1::QuonConst{S1}, p2::QuonConst{S2}) where {S1, S2} = (p1 === p2) || (change_direction(p1) === p2)
# Base.:(==)(p1::QuonParam, ::typeof(ParaHalfPi)) = 

const eq_cd = [
    (ParaHalfPi, PerpNegHalfPi),
    (PerpHalfPi, ParaNegHalfPi),
    (ParaNegHalfPi, PerpHalfPi),
    (PerpNegHalfPi, ParaHalfPi),
    (ParaPi, PerpInfPi),
    (PerpPi, ParaInfPi),
    (ParaZero, PerpInfZero),
    (PerpZero, ParaInfZero),
    (ParaInfZero, PerpZero),
    (PerpInfZero, ParaZero),
    (ParaInfPi, PerpPi),
    (PerpInfPi, ParaPi),
]
for (p1, p2) in eq_cd
    @eval change_direction(::typeof($p1)) = $p2
end

is_singular_change_direction(p::QuonConst{S}) where S = false
is_singular_yang_baxter(p1::QuonConst{S1}, p2::QuonConst{S2}, p3::QuonConst{S3}) where {S1, S2, S3} = false
for p in (ParaHalfPi, ParaNegHalfPi, ParaPi, ParaZero, ParaInfZero, ParaInfPi)
    @eval is_parallel(::typeof($p)) = true
end
for p in (PerpHalfPi, PerpNegHalfPi, PerpPi, PerpZero, PerpInfZero, PerpInfPi)
    @eval is_parallel(::typeof($p)) = false
end

is_para_pi(p::QuonConst) = (p == ParaPi)
is_perp_pi(p::QuonConst) = (p == ParaPi)
is_para_zero(p::QuonConst) = (p == ParaZero)
is_perp_zero(p::QuonConst) = (p == PerpZero)
is_para_half_pi(p::QuonConst) = (p == ParaHalfPi)
is_perp_half_pi(p::QuonConst) = (p == PerpHalfPi)

# function yang_baxter_param(p1::QuonConst, p2::QuonConst, p3::QuonConst)
#     # triangle => star
#     p1.isparallel && (p1 = change_direction(p1))
#     !p2.isparallel && (p2 = change_direction(p2))
#     p3.isparallel && (p3 = change_direction(p3))
    
#     q1_param, q2_param, q3_param = yang_baxter_param(p1.param, p2.param, p3.param)
#     q1, q2, q3 = QuonConst(q1_param, true), QuonConst(q2_param, false), QuonConst(q3_param, true)
#     return q1, q2, q3
# end
# function yang_baxter_param_inv(p1::QuonConst, p2::QuonConst, p3::QuonConst)
#     # star => triangle
#     !p1.isparallel && (p1 = change_direction(p1))
#     p2.isparallel && (p2 = change_direction(p2))
#     !p3.isparallel && (p3 = change_direction(p3))

#     q1_param, q2_param, q3_param = yang_baxter_param_inv(p1.param, p2.param, p3.param)
#     return QuonConst(q1_param, false), QuonConst(q2_param, true), QuonConst(q3_param, false)
# end
