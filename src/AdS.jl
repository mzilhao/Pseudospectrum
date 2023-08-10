module AdS

using LinearAlgebra
using ..Pseudospectrum

@inline r_of_x(x) = tan(pi * x / 2)
@inline V_of_r(r, ℓ::Int) = (r^2 + 1) * (2 + ℓ*(ℓ+1) / r^2)
@inline function V_of_x(x, ℓ::Int)
    rr = r_of_x(x)
    V_of_r(rr, ℓ)
end

function build_operator_AdS4_sph(T::Type, N::Int, ℓ::Int)
    x, D, D2_ = cheb(T(0), T(1), N)
    M = N - 1

    Id = Matrix(I, M, M)
    zero_mat = zeros(Int8, M, M)

    # by sampling only the points in x[2:end-1], we're removing the points x = 0
    # and x = 1. that's also why we're considering M = N - 1 above (remember
    # that length(x) = N+1).
    V_Id = V_of_x.(x[2:end-1], ℓ) .* Id
    D2 = D2_[2:end-1, 2:end-1]

    C = 4 / pi^2 * D2 - V_Id

    #= build 2M x 2M matrix

     ( 0  |  4/π² ∂_xx - V Id )
     ( Id |  0                )

    =#
    L = [[zero_mat C]; [Id zero_mat]]

    im .* L
end

struct AdS4_sph{S1,S2,S3} <: AbstractOperator
    L :: S1
    G :: S2
    A :: S3
end

function AdS4_sph(T::Type, N::Int, ℓ::Int)
    L  = build_operator_AdS4_sph(T, N, ℓ)
    G_ = build_Gram_matrix_AdS4_sph(T, N, ℓ)

    #=
    G should be Hermitian by construction (TODO: add tests for this!), but due
    to round-off errors it won't be exactly. so force it to be thus, as
    mentioned in LinearAlgebra/src/cholesky.jl:361 (or, equivalently, the
    documentation for the cholesky method:

      If you have a matrix A that is slightly non-Hermitian due to roundoff
      errors in its construction, wrap it in Hermitian(A) before passing it to
      cholesky in order to treat it as perfectly Hermitian.
    =#
    G = Hermitian(G_)
    F = cholesky(G)
    A = F.U * L * inv(F.U)

    AdS4_sph(L,G,A)
end

AdS4_sph(N::Int, ℓ::Int) = AdS4_sph(Float64, N, ℓ)


function build_Gram_matrix_AdS4_sph(T::Type, N::Int, ℓ::Int)
    M = 2*N

    xN, D, _ = cheb(T(0), T(1), N)
    xM, _    = cheb(T(0), T(1), M)

    zero_mat = zeros(Int8, N-1, N-1)

    # interpolation matrix for the (M+1)-point grid
    II = interp_matrix(M,N)

    # integration weights in the (M+1)-point grid
    w = clencurt(T(0), T(1), M)

    CM1 = Diagonal(w)
    CM2 = 4/T(pi)^2 * Diagonal(w)

    CMV = V_of_x.(xM, ℓ) .* Diagonal(w)
    # remove singular points, which add to zero anyway
    CMV[1,:]   .= 0
    CMV[end,:] .= 0

    CN1 = II' * CM1 * II
    CN2 = II' * CM2 * II
    CNV = II' * CMV * II

    GE1_ = T(pi)/4 .* CN1
    GE2_ = T(pi)/4 .* (CNV .+ D' * CN2 * D)

    # removing the points x = 0 and x = 1, since these add to zero
    GE1 = GE1_[2:end-1, 2:end-1]
    GE2 = GE2_[2:end-1, 2:end-1]

    #= build 2(N-1) x 2(N-1) matrix

     ( GE1 | 0   )
     ( 0   | GE2 )

    =#
    GE = [[GE1 zero_mat];
          [zero_mat GE2]]

    GE
end


end
