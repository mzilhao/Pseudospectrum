module AdS5_IEF

using LinearAlgebra
using ..Pseudospectrum

function build_operators(T::Type, N::Int, k2::Number)
    z, D, D2 = cheb(T(0), T(1), N)

    N  = length(z)
    Id = Matrix(I, N, N)

    z2 = z .* z
    z3 = z .* z2
    z4 = z2 .* z2

    iL1  = @. z * (1 - z4) * D2 + (5 - 9z4) * D - (k2 * z + 16z3) * Id

    L2  = @. 2 * z * D + 5*Id

    im * iL1, L2
end

function build_Gram_matrix(T::Type, N::Int, k2::Number)
    M = 2*N

    zN, D, _ = cheb(T(0), T(1), N)
    zM, _    = cheb(T(0), T(1), M)

    # interpolation matrix for the (M+1)-point grid
    II = interp_matrix(M,N)

    # integration weights in the (M+1)-point grid
    w = clencurt(T(0), T(1), M)
    W = Diagonal(w)

    CM1 = @. W * zM^3 * (k2 * zM^2 + 10 - 4*zM^4)
    CM2 = @. W * zM^5 * (1 - zM^4)
    CM3 = @. W * 5//2 * zM^4 * (1 - zM^4)

    CN1 = II' * CM1 * II
    CN2 = II' * CM2 * II
    CN3 = II' * CM3 * II

    GE = (CN1 .+ D' * CN2 * D .+ CN3 * D .+ D' * CN3) / 2

    GE
end


struct AdS5_planar_IEF{S1,S2,S3}
    A :: S1
    B :: S2
    G :: S3
end

function AdS5_planar_IEF(T::Type, N::Int, k2::Number)
    L1, L2  = build_operators(T, N, k2)
    G_ = build_Gram_matrix(T, N, k2)

    #=
    G should be Hermitian by construction, but due to round-off errors it won't
    be exactly. so force it to be thus, as mentioned in
    LinearAlgebra/src/cholesky.jl:361 (or, equivalently, the documentation for
    the cholesky method):

      If you have a matrix A that is slightly non-Hermitian due to roundoff
      errors in its construction, wrap it in Hermitian(A) before passing it to
      cholesky in order to treat it as perfectly Hermitian.
    =#
    G = Hermitian(G_)
    F = cholesky(G)
    A = F.U * L1 * inv(F.U)
    B = F.U * L2 * inv(F.U)

    AdS5_planar_IEF(A,B,G)
end
AdS5_planar_IEF(N::Int, k2::Number) = AdS5_planar_IEF(Float64, N, k2)

end
