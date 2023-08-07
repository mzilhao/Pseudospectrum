
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

@testset "Clenshaw-Curtis weights:" begin
    w10 = PS.clencurt(10)
    w10_true = [0.010101, 0.094579, 0.185635, 0.253588, 0.299213, 0.313766,
                0.299213, 0.253588, 0.185635, 0.094579, 0.010101]
    @test isapprox(w10, w10_true, atol=1e-6)

    w9 = PS.clencurt(9)
    w9_true = [0.012346, 0.116567, 0.225284, 0.301940, 0.343863, 0.343863,
               0.301940, 0.225284, 0.116567, 0.012346]
    @test isapprox(w9, w9_true, atol=1e-5)
end

@testset "Clenshaw-Curtis integration:" begin
    N = 2
    x, = PS.cheb(N)

    f = -x.^2 .+ 1
    w = PS.clencurt(N)

    int = w' * f
    @test int ≈ 4/3


    N = 3
    x, = PS.cheb(0.0, 1.0, N)

    f = 4 .* x.^3
    w = PS.clencurt(0.0, 1.0, N)

    int = w' * f
    @test int ≈ 1.0


    N = 12
    x, = PS.cheb(-pi/2, pi/2, N)

    f = cos.(x)
    w = PS.clencurt(-pi/2, pi/2, N)

    int = w' * f
    @test int ≈ 2.0
end
