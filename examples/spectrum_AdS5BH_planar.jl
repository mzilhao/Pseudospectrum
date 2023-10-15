using Pseudospectrum
using LinearAlgebra

ksq = 0
N = 128

model = AdS5BH_planar(k2)
Op = Operators(N, model)

F = eigen(Op.A, Op.B; sortby = x -> real(x)^2 + imag(x)^2)

# follow the rule of thumb by Boyd and consider only the N/2 first eigenmodes.
# also, take the conjugate to remove the minus sign in all the imaginary parts
Nmax  = Int(N/2)
omega = conj(F.values[1:Nmax])


# using the Schur decomposition:

# F = schur(Op.A, Op.B)

# λ = F.α ./ F.β
# conj!(λ)
# sort!(λ, by=imag)


using Plots
gr()

scatter(real(omega), imag(omega))
