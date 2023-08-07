
@testset "cheb:" begin
    x, D = PS.cheb(2)
    @test D ≈ [[-1.5   2.0  -0.5];
               [-0.5   0.0   0.5];
               [ 0.5  -2.0   1.5]]

    x, D = PS.cheb(Float32, 2)
    @test D ≈ [[-1.5   2.0  -0.5];
               [-0.5   0.0   0.5];
               [ 0.5  -2.0   1.5]]
end
