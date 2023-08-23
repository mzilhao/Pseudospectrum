using Pseudospectrum
using LinearAlgebra

function build_operators(z::AbstractVector, D::AbstractMatrix, D2::AbstractMatrix,
                         k2::Number)
    N  = length(z)
    Id = Matrix(I, N, N)

    z2 = z .* z
    z4 = z2 .* z2
    z5 = z .* z4

    L1 = @. -(1 - z4) * z2 * D2 +
        4*z5 * D +
        ( z2 * k2 + 3//4 * (1 - z4) + 3 * (1 + z4) ) * Id

    L2 = @. 2 * z2 * im * D

    L1, L2
end

ksq = 0
N = 128

T = Float64

z, D, D2 = cheb(zero(T), one(T), N)
L1, L2 = build_operators(z, D, D2, ksq)
F = eigen(L1, L2; sortby = x -> -imag(x))


# follow the rule of thumb by Boyd and consider only the N/2 first eigenmodes.
# also, take the conjugate to remove the minus sign in all the imaginary parts
Nmax  = Int(N/2)
omega = conj(F.values[1:Nmax])


# using the Schur decomposition:

# F = schur(L1, L2)

# λ = F.α ./ F.β
# conj!(λ)
# sort!(λ, by=imag)


using Plots
gr()

scatter(real(omega), imag(omega))
