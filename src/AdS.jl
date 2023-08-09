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

end
