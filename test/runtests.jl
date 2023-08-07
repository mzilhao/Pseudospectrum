using Pseudospectrum
using Test

const PS = Pseudospectrum

@testset "Pseudospectrum.jl" begin
    @time @testset "utils tests:" begin include("utils_tests.jl") end
end
