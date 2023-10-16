
@testset "AdS5BH (IEF) spectrum test:" begin
    N = 64

    k2 = 0
    model = AdS5BH_planar(k2)
    Op = Operators(N, model)

    F = eigen(Op.A, Op.B; sortby = x -> real(x)^2 + imag(x)^2)

    # let's compare the first 2 values, which should be ω ≈ ±3.119451591 - i 2.746675744
    ω1 = F.values[1]
    ω2 = F.values[2]

    Reω_target = 3.119451591
    Imω_target = -2.746675744

    @test real(ω1) ≈ Reω_target || real(ω1) ≈ -Reω_target
    @test imag(ω1) ≈ Imω_target

    @test real(ω1) ≈ -real(ω2)
    @test imag(ω1) ≈ imag(ω2)

    # now the following 2 values, which should be ω ≈ ±5.169521 - i 4.763570
    ω3 = F.values[3]
    ω4 = F.values[4]

    Reω_target = 5.169521
    Imω_target = -4.763570

    @test isapprox(real(ω3), Reω_target, atol=1e-6)  || isapprox(real(ω3), -Reω_target, atol=1e-6)
    @test isapprox(imag(ω3), Imω_target, atol=1e-6)

    @test real(ω3) ≈ -real(ω4)
    @test imag(ω3) ≈ imag(ω4)


    # and now for the following 2, which should be ω ≈ ±7.188 - i 6.769
    ω5 = F.values[5]
    ω6 = F.values[6]

    Reω_target = 7.188
    Imω_target = -6.769

    @test isapprox(real(ω5), Reω_target, atol=1e-3)  || isapprox(real(ω5), -Reω_target, atol=1e-3)
    @test isapprox(imag(ω6), Imω_target, atol=1e-3)

    @test isapprox(real(ω5), -real(ω6), atol=1e-5)
    @test isapprox(imag(ω5), imag(ω6), atol=1e-5)
end

@testset "test if Gram matrix Hermitian:" begin
    N = 12

    k2 = 0
    model = AdS5BH_planar(k2)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    k2 = 1
    model = AdS5BH_planar(k2)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    k2 = 4
    model = AdS5BH_planar(k2)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    k2 = 8
    model = AdS5BH_planar(k2)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'
end

@testset "AdS5BH (IEF) pseudospectrum test:" begin
    N = 64

    k2 = 0
    model = AdS5BH_planar(k2)
    Op = Operators(N, model)

    # QNMs
    z = [3.119451591 - 2.746675744im, -3.119451591 - 2.746675744im, 5.169521 - 4.763570im, -5.169521 - 4.763570im,
         7.188 - 6.7696im, -7.188 - 6.7696im]

    sigmin = [basic_svd(real(zi), imag(zi), Op.A, Op.B) for zi in z]
    @test all(sigmin .< 1e-10)
end
