using Pseudospectrum
using LinearAlgebra

using Plots
using LaTeXStrings


ℓ = 4
N = 32

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

sigmin = basic_svd(x, y, Op.A)

# spectra
F = eigen(Op.L; sortby = x -> abs(real(x)))
omega = F.values

#Nmax = Int(floor(N/2))
Nmax = 10
ω   = omega[1:Nmax]
Reω = real(ω)
Imω = imag(ω)


plotlyjs()
#pyplot()

#gr()



contourf(x, y, sigmin', color=:plasma, aspect_ratio=:equal, levels=12,
         title=L"AdS $\ell=2$")

scatter!(Reω, Imω, shape=:x, color="red", markersize=3, label="Normal Modes")

xlabel!(L"$\Re(\omega)$")
ylabel!(L"$\Im(\omega)$")

#savefig("AdS4_PS.pdf")
