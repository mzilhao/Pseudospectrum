module AdS_IEF

using LinearAlgebra
using ..Pseudospectrum
using ..AdS: V_of_x

function build_operator_AdS4_sph(T::Type, N::Int, ℓ::Int)
    x, D_, _ = cheb(T(0), T(1), N)
    M = N - 1

    Id = Matrix(I, M, M)
    zero_mat = zeros(Int8, M, M)

    ## TODO: the procedure below cannot accommodate the case ℓ=0

    # by sampling only the points in x[2:end-1], we're removing the points x = 0
    # and x = 1. that's also why we're considering M = N - 1 above (remember
    # that length(x) = N+1).
    V_Id = V_of_x.(x[2:end-1], ℓ) .* Id
    D = D_[2:end-1, 2:end-1] ./ T(pi)

    #= build 2M x 2M matrix

     ( -2/π ∂_x  | Id       )
     ( V/2 Id    | -1/π ∂_x )

    =#
    L = [[-2*D     Id];
         [ V_Id/2  -D]]

    im * L
end

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
    CMV = V_of_x.(xM, ℓ) .* Diagonal(w)
    # make the singular points zero, since they add to zero anyway due to the
    # boundary conditions
    CMV[1,:]   .= 0
    CMV[end,:] .= 0

    CN1 = II' * CM1 * II
    CNV = II' * CMV * II

    GE1_ = T(pi)/4 .* CNV
    GE2_ = T(pi)/4 .* CN1

    ## TODO: the procedure below cannot accommodate the case ℓ=0

    # remove the rows and columns multiplying with points x = 0 and x = 1,
    # since these add to zero
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


struct AdS4_sph_IEF{S1,S2,S3,S4}
    L :: S1
    G :: S2
    A :: S3
    B :: S4
end

function AdS4_sph_IEF(T::Type, N::Int, ℓ::Int)
    @assert ℓ>=0
    @assert ℓ != 0 "Case ℓ=0 is not yet correctly implemented"
    L  = build_operator_AdS4_sph(T, N, ℓ)
    G_ = build_Gram_matrix_AdS4_sph(T, N, ℓ)

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
    A = F.U * L * inv(F.U)

    #= build 2(N-1) x 2(N-1) matrix

     ( 0 | 0  )
     ( 0 | Id )

    =#
    zero_mat = zeros(Int8, N-1, N-1)
    Id = Matrix(I, N-1, N-1)
    B = [[zero_mat zero_mat];
         [zero_mat Id]]

    AdS4_sph_IEF(L,G,A,B)
end
AdS4_sph_IEF(N::Int, ℓ::Int) = AdS4_sph_IEF(Float64, N, ℓ)

end
