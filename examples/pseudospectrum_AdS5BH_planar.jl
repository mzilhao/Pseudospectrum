using Pseudospectrum
using LinearAlgebra

using Plots
using LaTeXStrings


k2 = 0
N  = 48

model = AdS5BH_planar(k2)
Op = Operators(N, model)

# spectra
F = eigen(Op.A, Op.B; sortby = x -> real(x)^2 + imag(x)^2)
omega = F.values


function plot_ps(x, y, Op)
    sigmin = basic_svd(x, y, Op.A, Op.B)

    p = contourf(x, y, log10.(sigmin'), color=:viridis, aspect_ratio=:equal,
             levels=16,
             title=L"AdS5 BH $k^2=0$")

    xlabel!(L"$\Re(\omega)$")
    ylabel!(L"$\Im(\omega)$")

    return p
end


xmin = -10.0
xmax =  10.0
ymin = -12.0
ymax =   0.0

Nx = 201
Ny = 201
x  = range(xmin, xmax, length=Nx)
y  = range(ymin, ymax, length=Ny)

plot_obj = plot_ps(x, y, Op)

# add QNMs
Nmax = 8
ω    = omega[1:Nmax]
Reω  = real(ω)
Imω  = imag(ω)

scatter!(Reω, Imω, shape=:x, color="red", markersize=3, label="QNM")

savefig("AdS5BH_PS.pdf")



# zoom in on first mode
xmin =  2.9
xmax =  3.5
ymin = -3.0
ymax = -2.4

Nx = 301
Ny = 301
x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

plot_obj = plot_ps(x, y, Op)

# add QNM
ω   = [omega[1]]
Reω = real(ω)
Imω = imag(ω)

scatter!(Reω, Imω, shape=:x, color="red", markersize=3, label="")

savefig("AdS5BH_PS_zoom.pdf")
