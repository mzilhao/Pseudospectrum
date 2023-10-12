using Pseudospectrum
using LinearAlgebra

using Plots
using LaTeXStrings


k2 = 0
N  = 48

Op = AdS5_planar_IEF(N, k2)


Nx = 201
Ny = 201

T = Float64

xmin = -T(10)
xmax =  T(10)

ymin = -T(12)
ymax =  T(0)

x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

sigmin = basic_svd(x, y, Op.A, Op.B)

# spectra
F = eigen(Op.A, Op.B; sortby = x -> real(x)^2 + imag(x)^2)
omega = F.values

#Nmax = Int(floor(N/2))
Nmax = 8
ω   = omega[1:Nmax]
Reω = real(ω)
Imω = imag(ω)


#plotlyjs()
#pyplot()

gr()



# contourf(x, y, sigmin', color=:plasma, aspect_ratio=:equal,
#          levels=12,
#          title=L"AdS5 $k^2=0$")

contourf(x, y, log10.(sigmin'), color=:viridis, aspect_ratio=:equal,
         levels=12,
         title=L"AdS5 $k^2=0$")

scatter!(Reω, Imω, shape=:x, color="red", markersize=3, label="QNM")

xlabel!(L"$\Re(\omega)$")
ylabel!(L"$\Im(\omega)$")

savefig("AdS5_IEF_PS.png")
