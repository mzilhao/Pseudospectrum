
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

@testset "Clenshaw-Curtis:" begin
    w10 = PS.clencurt(10)
    w10_true = [0.010101, 0.094579, 0.185635, 0.253588, 0.299213, 0.313766,
                0.299213, 0.253588, 0.185635, 0.094579, 0.010101]
    @test isapprox(w10, w10_true, atol=1e-6)

    w9 = PS.clencurt(9)
    w9_true = [0.012346, 0.116567, 0.225284, 0.301940, 0.343863, 0.343863,
               0.301940, 0.225284, 0.116567, 0.012346]
    @test isapprox(w9, w9_true, atol=1e-5)
end
