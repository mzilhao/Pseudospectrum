using LinearAlgebra

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

@testset "Interpolation:" begin
    # for rectangular cases (equal grids), the interpolation matrix should
    # reduce to the identity matrix.
    II = PS.interp_matrix(5,5)
    @test II ≈ Matrix(I, 6, 6)

    II = PS.interp_matrix(32,32)
    @test II ≈ Matrix(I, 33, 33)


    # upsample, from a grid with size 5+1 to one with size 7+1
    N = 5
    x, = PS.cheb(N)
    # "old" vector
    fN = 1 .- x.^2 .+ x.^3

    # new grid
    M = 7
    xnew, = PS.cheb(M)

    # interpolation matrix
    II = PS.interp_matrix(M, N)

    # new vector
    fM = II * fN
    @test fM ≈ 1 .- xnew.^2 .+ xnew.^3


    # downsample, from a grid with size 12+1 to one with size 5+1
    N = 12
    x, = PS.cheb(N)
    # "old" vector
    fN = x .- x.^4 .+ 42

    # new grid
    M = 5
    xnew, = PS.cheb(M)

    # interpolation matrix
    II = PS.interp_matrix(M, N)

    # new vector
    fM = II * fN
    @test fM ≈ xnew .- xnew.^4 .+ 42


    # more examples

    N = 24
    x, = PS.cheb(N)
    # "old" vector
    fN = sin.(2*pi .* x) .+ 0.1 * cos.(pi .* x)

    # new grid
    M = 32
    xnew, = PS.cheb(M)

    # interpolation matrix
    II = PS.interp_matrix(M, N)

    # new vector
    fM = II * fN
    @test fM ≈ sin.(2*pi .* xnew) .+ 0.1 * cos.(pi .* xnew)

end
