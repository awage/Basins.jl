using Revise
using Plots
using DynamicalSystems
using Basins


# Simplest Chaotic Flows withInvolutional SymmetriesJ. C. Sprott IJBC 2013
@inline @inbounds function dl_lorenz(u, p, t)
    R = p[1];
	x = u[1]; y = u[2]; z = u[3];
    dx = y - x
    dy = -x*z
	dz = x*y - R
    return SVector{3}(dx, dy, dz)
end


R=4.7;
p= [R]
ds = ContinuousDynamicalSystem(dl_lorenz, rand(3), p)
xg=range(-10.,10.,length=400)
yg=range(-10.,10.,length=400)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-12, atol = 1e-12), reltol=1e-10, abstol=1e-10)
@time bsn = Basins.basins_map2D(xg, yg, pmap; Ncheck=10)
#integ = integrator(ds, alg=Tsit5(), reltol=1e-9, abstol=1e-9, save_everystep=false)
#@time bsn = Basins.basin_general_ds(xg, yg, integ; dt=2., idxs=1:2)

plot(xg,yg,bsn.basin',seriestype=:heatmap)
