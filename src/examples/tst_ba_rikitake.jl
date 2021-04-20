using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
integ=integrator(ds)
xg=range(-6.,6.,length=300)
yg=range(-6.,6.,length=300)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9)
@time bsn = basin_poincare_map(xg, yg, pmap)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
