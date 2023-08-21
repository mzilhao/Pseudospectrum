using Pseudospectrum
using LinearAlgebra

ℓ = 2
N = 64
#N = 12

# Op = AdS4_sph(N, ℓ)

# F = eigen(Op.L; sortby = x -> abs(real(x)))
# omega = F.values


Op_IEF = AdS4_sph_IEF(N, ℓ)

F = eigen(Op_IEF.L, Op_IEF.B; sortby = x -> abs(real(x)))
omega = F.values

