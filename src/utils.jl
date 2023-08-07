
"""
    cheb(xmin, xmax, N)

Return a Chebyshev-Lobatto grid, together with its first and second
differentiation matrices.

# Arguments
* `xmin::Real`: rightmost grid point
* `xmax::Real`: leftmost grid point
* `N::Integer`: total number of grid points
"""
function cheb(xmin::T, xmax::T, N::Integer) where {T<:Real}
    x, D, D2 = cheb(N)
    x = (xmax + xmin .+ (xmax - xmin) * x) / 2
    D  ./= (xmax - xmin) / 2
    D2 ./= (xmax - xmin)^2 / 4
    x, D, D2
end


# taken from Trefethen (2000), "Spectral Methods in MatLab"

"When omitting grid limits, default to [-1,1] interval"
function cheb(N::Integer)
    @assert(N > 0)
    x = -cos.(pi*(0:N)/N)

    c  = [2; ones(N-1, 1); 2] .* (-1).^(0:N)
    X  = repeat(x, 1, N+1)
    dX = X - X'

    D = (c * (1 ./ c)') ./ (dX + Matrix(I, N+1, N+1))  # off-diagonal entries
    D = D - diagm(0 => sum(D', dims=1)[:])             # diagonal entries

    x, D, D*D
end
