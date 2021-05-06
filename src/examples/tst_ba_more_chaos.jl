using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins
using ChaosTools


ds = Systems.more_chaos_example()
integ=integrator(ds)
xg=range(-10.,10.,length=100)
yg=range(-10.,10.,length=100)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9)
@time bsn = Basins.basins_map2D(xg, yg, pmap)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
