using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

### First version using high-level function.
b=0.1665
ds = Systems.thomas_cyclical(b = b)
integ=integrator(ds, reltol=1e-9)
xg=range(-6.,6.,length=400)
yg=range(-6.,6.,length=400)
@time basin = draw_basin(xg, yg, integ; plane=(3, 0.), idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8))
plot(xg,yg,basin',seriestype=:heatmap)
