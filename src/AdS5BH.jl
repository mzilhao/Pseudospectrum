
struct AdS5BH_planar{T<:Number} <: Model
    k2 :: T
end

function build_operators(T::Type, N::Int, model::AdS5BH_planar)
    k2 = model.k2

    z, D, D2 = cheb(T(0), T(1), N)

    N  = length(z)
    Id = Matrix(I, N, N)

    z2 = z .* z
    z3 = z .* z2
    z4 = z2 .* z2

    iL1 = @. z * (1 - z4) * D2 + (5 - 9z4) * D - (k2 * z + 16z3) * Id
    L2  = @. 2 * z * D + 5*Id

    im * iL1, L2
end

function build_Gram_matrix(T::Type, N::Int, model::AdS5BH_planar)
    k2 = model.k2

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


    # note: the following is supposed to be a better way to perform this
    # computation, since it more "correctly" accounts for the interpolation
    # matrices, but in practice results with that method look more "noisy" for
    # higher overtones. so we leave it commented out, and for further
    # investigation in the future if needed.

    # zN, _    = cheb(T(0), T(1), N)
    # zM, D, _ = cheb(T(0), T(1), M)

    # # interpolation matrix for the (M+1)-point grid
    # II = interp_matrix(M,N)

    # # integration weights in the (M+1)-point grid
    # w = clencurt(T(0), T(1), M)
    # W = Diagonal(w)

    # CM1 = @. W * zM^3 * (k2 * zM^2 + 10 - 4*zM^4)
    # CM2 = @. W * zM^5 * (1 - zM^4)
    # CM3 = @. W * 5//2 * zM^4 * (1 - zM^4)

    # GE1 = CM1
    # GE2 = D' * CM2 * D
    # GE3 = CM3 * D

    # CN1 = II' * GE1 * II
    # CN2 = II' * GE2 * II
    # CN3 = II' * GE3 * II

    # GE = (CN1 .+ CN2 .+ CN3 .+ CN3') / 2


    GE
end
