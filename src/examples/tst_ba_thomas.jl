using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

### First version using high-level function.
b=0.1665
#b=0.22
ds = Systems.thomas_cyclical(b = b)
#integ=integrator(ds, reltol=1e-9)
xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9)

@time bsn = Basins.basins_map2D(xg, yg, pmap)

plot(xg,yg,bsn.basin',seriestype=:heatmap)
