
"""
    basic_svd(x, y, A)

Compute pseudospectra using a basic SVD decomposition. Based on page 371 (section 39) from
Trefethen, Embree, "Spectra and Pseudospectra - The Behavior of Nonnormal Matrices and Operators"

"""
@inline function basic_svd(x::Number, y::Number, A::AbstractMatrix)
    λ = x + im * y
    M = λ * I - A
    svd_M = svd(M)
    # the singular values in S are sorted in descending order, so the last
    # element is the smaller one, which is what we want
    svd_M.S[end]
end

function basic_svd(x::AbstractArray, y::AbstractArray, A::AbstractMatrix)
    T = eltype(x)
    Nx = length(x)
    Ny = length(y)
    sigmin = zeros(T,Nx,Ny)
    basic_svd!(sigmin, x, y, A)
end

"""
    basic_svd!(sigmin, x, y, A)
"""
function basic_svd!(sigmin, x::AbstractArray, y::AbstractArray, A::AbstractMatrix)
    @assert size(sigmin) == length(x), length(y)

    @inbounds for i in eachindex(x)
        @inbounds for j in eachindex(y)
            sigmin[i,j] = basic_svd(x[i], y[i], A)
        end
    end

    sigmin
end
