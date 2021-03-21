using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins
using ChaosTools

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
iter_f!,integ = poincaremap(ds, (3, 0.), 20., direction=+1, idxs=[1,2], rootkw = (xrtol = 1e-8, atol = 1e-8),reltol=1e-9)
reinit_f! = (integ,y) -> reinit!(integ,[y...,0.], t0=0.)

xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)

@time basin=draw_basin(xg, yg, integ, iter_f!, reinit_f!)

plot(xg,yg,basin',seriestype=:heatmap)

reinit_f!(integ, [5.1, -2.5])
u=Dataset([iter_f!(integ) for k=1:1000])
plot(u[600:1000,1], u[600:1000,2], seriestype=:scatter)
