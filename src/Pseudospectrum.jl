module Pseudospectrum

using LinearAlgebra

abstract type AbstractOperator end

include("utils.jl")
include("AdS.jl")

export AbstractOperator
export AdS4_sph

end
