module Pseudospectrum

using LinearAlgebra

include("utils.jl")
export cheb, clencurt, interp_matrix

include("compute.jl")
export basic_svd, basic_svd!

include("AdS.jl")
using .AdS: AdS4_sph
export AdS4_sph

include("AdS_IEF.jl")
using .AdS_IEF: AdS4_sph_IEF
export AdS4_sph_IEF

end
