using Pseudospectrum
using LinearAlgebra

using Plots
using LaTeXStrings


ℓ = 2
N = 24

Op = AdS4_sph_IEF(N, ℓ)


Nx = 301
Ny = 31

T = Float64

xmin = -T(15)
xmax =  T(15)

ymin = -T(5)
ymax =  T(5)

x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

sigmin = basic_svd(x, y, Op.A, Op.B)

# spectra
F = eigen(Op.A, Op.B; sortby = x -> abs(real(x)))
omega = F.values

#Nmax = Int(floor(N/2))
Nmax = 9
ω   = omega[2:Nmax]
Reω = real(ω)
Imω = imag(ω)


#plotlyjs()
#pyplot()

gr()



contourf(x, y, sigmin', color=:plasma, aspect_ratio=:equal, levels=12,
         title=L"AdS $\ell=2$")

scatter!(Reω, Imω, shape=:x, color="red", markersize=3, label="Normal Modes")

xlabel!(L"$\Re(\omega)$")
ylabel!(L"$\Im(\omega)$")

#savefig("AdS4_IEF_PS.pdf")
