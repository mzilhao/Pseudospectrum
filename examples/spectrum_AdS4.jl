using Pseudospectrum
using LinearAlgebra

ll = 2
N = 64

Op = AdS4_sph(N, ll)


F = eigen(Op.L; sortby = x -> abs(real(x)))
omega = F.values

