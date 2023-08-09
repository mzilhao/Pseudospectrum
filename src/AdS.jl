module AdS

using LinearAlgebra
using ..Pseudospectrum

@inline r_of_x(x) = tan(pi * x / 2)
@inline V_of_r(r, ll::Int) = (r^2 + 1) * (2 + ll*(ll+1) / r^2)
@inline function V_of_x(x, ll::Int)
    rr = r_of_x(x)
    V_of_r(rr, ll)
end

function build_operator_AdS4_sph(x::AbstractVector, D::AbstractMatrix,
                                 D2_::AbstractMatrix, ll::Int)
    N  = length(x) - 2
    Id = Matrix(I, N, N)
    zero_mat = zeros(Int8, N, N)

    # by sampling only the points in x[2:end-1], we're removing the points x = 0
    # and x = 1. that's also why we're considering N = length(x) - 2 above.

    V_Id = V_of_x.(x[2:end-1], ll) .* Id
    D2 = D2_[2:end-1, 2:end-1]

    C = 4 / pi^2 * D2 - V_Id

    #= build 2N matrix

     ( 0  |  4/π² ∂_xx - V Id )
     ( Id |  0                )

    =#
    L = [[zero_mat C]; [Id zero_mat]]

    im .* L
end

struct AdS4_sph{S} <: AbstractOperator
    L :: S
end

function AdS4_sph(T::Type, N::Int, ll::Int)
    x, D, D2 = cheb(T(0), T(1), N)
    L  = build_operator_AdS4_sph(x, D, D2, ll)

    AdS4_sph(L)
end

AdS4_sph(N::Int, ll::Int) = AdS4_sph(Float64, N, ll)


# TODO: mudar ll ou II
function build_Gram_matrix_AdS4_sph(T::Type, N::Int, ll::Int)
    M = 2*N

    xN, D, _ = cheb(T(0), T(1), N)
    xM, _    = cheb(T(0), T(1), M)

    zero_mat = zeros(Int8, N+1, N+1)

    # interpolation matrix for the (M+1)-point grid
    II = interp_matrix(M,N)

    # integration weights in the (M+1)-point grid
    w = clencurt(T(0), T(1), M)

    CM1 = Diagonal(w)
    CM2 = 4/T(pi)^2 * Diagonal(w)

    CMV = V_of_x.(xM, ll) .* Diagonal(w)
    # remove singular points
    CMV[1,:]   .= 0
    CMV[end,:] .= 0

    CN1 = II' * CM1 * II
    CN2 = II' * CM2 * II
    CNV = II' * CMV * II

    GE1 = T(pi)/4 .* CN1
    GE2 = T(pi)/4 .* (CNV .+ D' * CN2 * D)

    #= build 2(N+1) matrix

     ( GE1 | 0   )
     ( 0   | GE2 )

    =#
    GE = [[GE1 zero_mat];
          [zero_mat GE2]]

    GE
end


end
