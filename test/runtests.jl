using Pseudospectrum
using Test

const PS = Pseudospectrum

@testset "Pseudospectrum.jl" begin
    @time @testset "utils tests:" begin include("utils_tests.jl") end
    @time @testset "AdS4 tests:"   begin include("AdS4_tests.jl") end
    @time @testset "AdS5BH (IEF) tests:"  begin include("AdS5BH_tests.jl") end
end
