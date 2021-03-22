using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins
using ChaosTools

# reinit_f!(integ, [5.1, -2.5])
# u=Dataset([iter_f!(integ) for k=1:1000])
# plot(u[600:1000,1], u[600:1000,2], seriestype=:scatter)

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
integ=integrator(ds)
xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)
@time basin = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2)
plot(xg,yg,basin',seriestype=:heatmap)
