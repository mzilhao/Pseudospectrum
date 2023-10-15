module Pseudospectrum

using LinearAlgebra

abstract type Model end

# The following need a method for each `Model`

@doc raw"""
    build_operators(T::Type, N::Int, model::Model)

Build operators L1 and L2
```math
L_2 \partial_t u = i L_1 u
```
"""
function build_operators end

# TODO: documentation
"""
    build_Gram_matrix(T::Type, N::Int, model::Model)

"""
function build_Gram_matrix end


struct Operators{S1,S2,S3}
    A :: S1
    B :: S2
    G :: S3
end

function Operators(T::Type, N::Int, model::Model)
    L1, L2  = build_operators(T, N, model)
    G_ = build_Gram_matrix(T, N, model)

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

    Operators(A,B,G)
end
Operators(N::Int, model::Model) = Operators(Float64, N, model)

export Model, Operators


include("utils.jl")
export cheb, clencurt, interp_matrix

include("compute.jl")
export basic_svd, basic_svd!


include("AdS.jl")
export AdS4_sph

include("AdS5_IEF.jl")
export AdS5_planar_IEF

end
