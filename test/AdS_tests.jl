using LinearAlgebra

@testset "AdS spectrum test:" begin
    N = 64

    ℓ = 0

    Op = PS.AdS4_sph(N, ℓ)
    F  = eigen(Op.L; sortby = x -> abs(real(x)))
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

    Op = PS.AdS4_sph(N, ℓ)
    F  = eigen(Op.L; sortby = x -> abs(real(x)))
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

    Op = PS.AdS4_sph(N, ℓ)
    F  = eigen(Op.L; sortby = x -> abs(real(x)))
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
