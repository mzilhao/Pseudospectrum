using Pseudospectrum
using LinearAlgebra

ℓ = 2
N = 64

model = AdS4_sph(ℓ)

Op = Operators(N, model)

F = eigen(Op.A ; sortby = x -> abs(real(x)))
omega = F.values
