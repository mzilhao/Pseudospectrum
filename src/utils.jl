
"""
    cheb(xmin, xmax, N)

Return a Chebyshev-Lobatto grid, together with its first and second
differentiation matrices.

# Arguments
* `xmin :: AbstractFloat`: rightmost grid point
* `xmax :: AbstractFloat`: leftmost grid point
* `N    :: Integer`: total number of grid points
"""
function cheb(xmin::T, xmax::T, N::Integer) where {T<:AbstractFloat}
    x, D, D2 = cheb(T,N)
    x = (xmax + xmin .+ (xmax - xmin) * x) / 2
    D  ./= (xmax - xmin) / 2
    D2 ./= (xmax - xmin)^2 / 4
    x, D, D2
end


# taken from Trefethen (2000), "Spectral Methods in MatLab"

"When omitting grid limits, default to [-1,1] interval"
function cheb(T::Type, N::Integer)
    @assert(N > 0)
    x = -cos.(T(pi)*(0:N)/N)

    c  = [2; ones(T, N-1, 1); 2] .* (-1).^(0:N)
    X  = repeat(x, 1, N+1)
    dX = X - X'

    D = (c * (1 ./ c)') ./ (dX + Matrix(I, N+1, N+1))  # off-diagonal entries
    D = D - diagm(0 => sum(D', dims=1)[:])             # diagonal entries

    x, D, D*D
end
cheb(N::Integer) = cheb(Float64, N)

Tn(n::Int, x) = cos.(n .* acos.(x))


"""
    clencurt([T=Float64,] N)

Weights for the Clenshaw-Curtis quadrature for a Chebyshev-Lobatto grid. Taken
from Trefethen (2000), "Spectral Methods in MatLab".

# Arguments
* `xmin :: AbstractFloat`: rightmost grid point
* `xmax :: AbstractFloat`: leftmost grid point
* `N    :: Integer`: total number of grid points
"""
function clencurt(xmin::T, xmax::T, N::Integer) where {T<:AbstractFloat}
    w = clencurt(T, N)
    w * (xmax - xmin) / 2
end

function clencurt(T::Type, N::Integer)
    theta = T(pi) * (N:-1:0)/N

    w = zeros(T, N+1)

    if N % 2 == 0
        w[1] = w[N+1] = T(1) / (N^2-1)
        for i in 2:N
            s = zero(T)
            for k in 1:N/2-1
                s += 2 * cos(2*k*theta[i]) / (4*k^2 - 1)
            end
            s += cos(N*theta[i]) / (N^2 - 1)
            w[i] = 1 - s
        end
    else
        w[1] = w[N+1] = T(1) / N^2
        for i in 2:N
            s = zero(T)
            for k in 1:(N-1)/2
                s += 2 * cos(2*k*theta[i]) / (4*k^2 - 1)
            end
            w[i] = 1 - s
        end
    end
    w[2:N] *= 2/N

    w
end
clencurt(N::Integer) = clencurt(Float64, N)


"""
    interp_matrix([T=Float64,] M, N)

Interpolation matrix between two Chebyshev-Lobatto grids of sizes M+1 and N+1.
Given an interpolant polynomial with degree N, this matrix maps it to one with
degree M.
"""
function interp_matrix(T::Type, M::Int, N::Int)
    x1, = cheb(T,M)
    x2, = cheb(T,N)

    II = zeros(T, M+1, N+1)

    for i1 in 0:M
        for i2 in 0:N
            (i2 == 0 || i2 == N) ? α = 2 : α = 1
            s = T(1 + (-1)^i1 * (-1)^i2)
            for j in 1:N-1
                s += 2 * Tn(j, x2[i2+1]) * Tn(j, x1[i1+1])
            end
            II[i1+1, i2+1] = s / (N * α)
        end
    end

    II
end
interp_matrix(M::Int, N::Int) = interp_matrix(Float64, M, N)
