using Pseudospectrum
using LinearAlgebra

using Plots

ℓ = 2
N = 64

Op = AdS4_sph(N, ℓ)


Nx = 201
Ny = 21

T = Float64

xmin = -T(15)
xmax =  T(15)

ymin = -T(5)
ymax =  T(5)

x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)


sigmin = zeros(T,Nx,Ny)
it = 0
for i in eachindex(x)
    for j in eachindex(y)
        global it += 1
        λ = x[i] + im * y[j]
        if it % 1000 == 0
            @show λ
        end
        M = λ * I - Op.A
        svd_M = svd(M)
        # the singular values in S are sorted in descending order, so the last
        # element is the smaller one, which is what we want
        sigmin[i,j] = svd_M.S[end]
    end
end


# spectra
F = eigen(Op.L; sortby = x -> abs(real(x)))
omega = F.values

#Nmax = Int(floor(N/2))
Nmax = 10
ω   = omega[1:Nmax]
Reω = real(ω)
Imω = imag(ω)


plotlyjs()

contourf(x, y, sigmin', color=:plasma, aspect_ratio=:equal, levels=24)

scatter!(Reω, Imω, shape=:x, color="red", markersize=5)
