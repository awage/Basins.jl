using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins
using ChaosTools


ds = Systems.more_chaos_example()
iter_f!,integ = poincaremap(ds, (3, 0), 20., direction=+1, idxs=[1,2])
reinit_f! = (integ,y) -> reinit!(integ,[y...,0.], t0=0)

xg=range(-10.,10.,length=100)
yg=range(-10.,10.,length=100)

@time basin=draw_basin(xg, yg, integ, iter_f!, reinit_f!)

plot(xg,yg,basin',seriestype=:heatmap)
