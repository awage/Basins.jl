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
basin = draw_basin(xg, yg, integ; plane=(3, 0.), idxs = 1:2)
plot(xg,yg,basin',seriestype=:heatmap)
