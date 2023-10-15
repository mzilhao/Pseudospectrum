module AdS4

@inline r_of_x(x) = tan(pi * x / 2)
@inline V_of_r(r, ℓ::Int) = (r^2 + 1) * (2 + ℓ*(ℓ+1) / r^2)
@inline function V_of_x(x, ℓ::Int)
    rr = r_of_x(x)
    V_of_r(rr, ℓ)
end

end

struct AdS4_sph <: Model
    ℓ :: Int
end


function build_operators(T::Type, N::Int, model::AdS4_sph)
    ℓ = model.ℓ
    V = AdS4.V_of_x

    x, D, D2_ = cheb(T(0), T(1), N)
    M = N - 1

    Id = Matrix(I, M, M)
    zero_mat = zeros(Int8, M, M)

    # by sampling only the points in x[2:end-1], we're removing the points x = 0
    # and x = 1. that's also why we're considering M = N - 1 above (remember
    # that length(x) = N+1).
    V_Id = V.(x[2:end-1], ℓ) .* Id
    D2 = D2_[2:end-1, 2:end-1]

    C = 4 / T(pi)^2 * D2 - V_Id

    #= build 2M x 2M matrix

     ( 0  |  4/π² ∂_xx - V Id )
     ( Id |  0                )

    =#
    L = [[zero_mat C]; [Id zero_mat]]

    im .* L, I
end

function build_Gram_matrix(T::Type, N::Int, model::AdS4_sph)
    ℓ = model.ℓ
    V = AdS4.V_of_x

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

    CMV = V.(xM, ℓ) .* Diagonal(w)
    # make the singular points zero, since they add to zero anyway due to the
    # boundary conditions
    CMV[1,:]   .= 0
    CMV[end,:] .= 0

    CN1 = II' * CM1 * II
    CN2 = II' * CM2 * II
    CNV = II' * CMV * II

    GE1_ = T(pi)/4 .* CN1
    GE2_ = T(pi)/4 .* (CNV .+ D' * CN2 * D)

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
