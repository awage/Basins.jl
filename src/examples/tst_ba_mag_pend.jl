using Plots
using DynamicalSystems


# https://cdn.ima.org.uk/wp/wp-content/uploads/2020/03/Chaos-in-the-Magnetic-Pendulum-from-MT-April-2020.pdf
#ds = mag_pendulum(γ=1, d=0.5, α=0.175, ω=1., N=4)

ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=4)
xg=range(-2, 2,length = 400)
yg=range(-2, 2,length = 400)
default_diffeq = (reltol = 1e-9,  alg = Vern9())
@time bsn,att = basins_of_attraction((xg, yg), ds; diffeq = default_diffeq)


# Uncertainty exponent for these parameter and grid
v, e, bd = basins_fractal_dimension(bsn)
α = 2 - bd

plot(xg, yg, bsn',seriestype = :heatmap)
