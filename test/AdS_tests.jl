using LinearAlgebra

@testset "AdS4 spectrum test:" begin
    N = 64

    ℓ = 0
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    F  = eigen(Op.A; sortby = x -> abs(real(x)))
    omega  = F.values

    # let's compare the N/2 first values, as the error increases for
    # higher-frequency ones
    Nmax = Int(floor(N/2))
    ω = omega[1:Nmax]
    idx = real(ω) .> 0
    Reω = real(ω[idx])

    @test all(abs.(imag(ω)) .< 1e-13)
    @test all(Reω .≈ [2*n + 3 + ℓ for n in 0:length(Reω)-1])


    ℓ = 2
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    F  = eigen(Op.A; sortby = x -> abs(real(x)))
    omega  = F.values

    # let's compare the N/2 first values, as the error increases for
    # higher-frequency ones
    Nmax = Int(floor(N/2))
    ω = omega[1:Nmax]
    idx = real(ω) .> 0
    Reω = real(ω[idx])

    @test all(abs.(imag(ω)) .< 1e-13)
    @test all(Reω .≈ [2*n + 3 + ℓ for n in 0:length(Reω)-1])


    ℓ = 5
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    F  = eigen(Op.A; sortby = x -> abs(real(x)))
    omega  = F.values

    # let's compare the N/2 first values, as the error increases for
    # higher-frequency ones
    Nmax = Int(floor(N/2))
    ω = omega[1:Nmax]
    idx = real(ω) .> 0
    Reω = real(ω[idx])

    @test all(abs.(imag(ω)) .< 1e-13)
    @test all(Reω .≈ [2*n + 3 + ℓ for n in 0:length(Reω)-1])
end

@testset "test if Gram matrix Hermitian:" begin
    N = 12

    ℓ = 0
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    ℓ = 1
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    ℓ = 2
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'

    ℓ = 8
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    G = PS.build_Gram_matrix(Float64, N, model)
    @test G ≈ G'
end

@testset "AdS4 pseudospectrum test:" begin
    N = 64

    ℓ = 0
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    # actual eigenvalues
    x = [2*n + 3 + ℓ for n in 0:10]

    sigmin = [basic_svd(xi, 0, Op.A) for xi in x]
    @test all(sigmin .< 1e-12)


    ℓ = 1
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    # actual eigenvalues
    x = [2*n + 3 + ℓ for n in 0:10]

    sigmin = [basic_svd(xi, 0, Op.A) for xi in x]
    @show sigmin
    @test all(sigmin .< 1e-11)


    ℓ = 8
    model = AdS4_sph(ℓ)
    Op = Operators(N, model)

    # actual eigenvalues
    x = [2*n + 3 + ℓ for n in 0:10]

    sigmin = [basic_svd(xi, 0, Op.A) for xi in x]
    @test all(sigmin .< 1e-11)
end
