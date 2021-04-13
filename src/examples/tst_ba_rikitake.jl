using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
integ=integrator(ds)
xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)
@time bsn = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
