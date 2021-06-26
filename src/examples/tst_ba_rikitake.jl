using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

μ=0.46
ds = Systems.rikitake(μ = μ, α = 1.0)
integ=integrator(ds)
xg=range(-3.,3.,length=300)
yg=range(-3.,3.,length=300)
zg=range(-1.,1.,length=30)
pmap = poincaremap(ds, (3, 0.), Tmax=1e4; idxs = 1:2,  rootkw = (xrtol = 1e-6, atol = 1e-6), reltol=1e-6, abstol=1e-6)
@time bsn = Basins.basins_map2D(xg, yg, pmap)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
