using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins


# http://dx.doi.org/10.1155/2013/256092
# Multiple Coexisting Attractors and Hysteresis in the Generalized Ueda Oscillator
@inline @inbounds function ueda(u, p, t)
    b = p[1]; F = p[2]; ω = p[3]; α = p[4]
	x=u[1]; y=u[2]; z=u[3];
    dx = y
    dy = -b*y-sign(x)*abs(x)^α+F*sin(z)
    dz = ω
    return SVector{3}(dx, dy, dz)
end

ω=1
α=4
p= [0.05, 7.5, ω, α]
ds = ContinuousDynamicalSystem(ueda, rand(3), p)
integ  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)

xres=100
yres=100

xg = range(-2.,2.,length=xres)
yg = range(-7.,7.,length=yres)

@time basin = Basins.basins_map2D(xg, yg, integ; T=2π/ω)

plot(xg,yg,basin', seriestype=:heatmap)
