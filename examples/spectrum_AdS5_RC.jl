using Pseudospectrum
using LinearAlgebra

function build_operators(z::AbstractVector, D::AbstractMatrix, D2::AbstractMatrix,
                         k2::Number)
    N  = length(z)
    Id = Matrix(I, N, N)

    z2 = z  .* z
    z3 = z  .* z2
    z4 = z2 .* z2


    L1 = @. (1 - z4) * D2 - 4 * z3 * D - ( k2 + 9//4 * z2 + 15//4 / z2 ) * Id
    L1 ./= 1 .+ z4

    L2 = @. -2 * z4 * D - 4 * z3 * Id
    L2 ./= 1 .+ z4

    # remove the points at z = 0
    M = N - 1

    L1_ = L1[2:end, 2:end]
    L2_ = L2[2:end, 2:end]
    Id_ = Matrix(I, M, M)
    zero_mat = zeros(Int8, M, M)

    #= build 2(N-1) x 2(N-1) matrix

     ( 0  | Id )
     ( L1 | L2 )

    =#
    
    iL = [[zero_mat Id_];
          [L1_      L2_]]

    iL
end

ksq = 0
N = 128

T = Float64

z, D, D2 = cheb(zero(T), one(T), N)


iL = build_operators(z, D, D2, ksq)

F = eigen(im * iL; sortby = x -> real(x)^2 + imag(x)^2)




# follow the rule of thumb by Boyd and consider only the N/2 first eigenmodes.
# also, take the conjugate to remove the minus sign in all the imaginary parts
Nmax  = Int(N/2)
omega = conj(F.values[1:Nmax])


# # using the Schur decomposition:

# # F = schur(L1, L2)

# # λ = F.α ./ F.β
# # conj!(λ)
# # sort!(λ, by=imag)


using Plots
gr()

scatter(real(omega), imag(omega))
