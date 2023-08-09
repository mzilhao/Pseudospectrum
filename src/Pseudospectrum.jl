module Pseudospectrum

using LinearAlgebra

abstract type AbstractOperator end
export AbstractOperator

include("utils.jl")
export cheb, clencurt, interp_matrix

include("AdS.jl")
using .AdS: AdS4_sph
export AdS4_sph

end
