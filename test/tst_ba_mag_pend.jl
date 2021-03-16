using Revise
using Plots
using Basins
using DifferentialEquations
using DynamicalSystems


ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)

iter_f! = (x) -> step!(x, 2*pi/ω, true)
reinit_f! = (integ,y) -> reinit!(integ,[y...,0.,0.], t0=0)

xg=range(-4,4,length=300)
yg=range(-4,4,length=300)


@time basin=draw_basin(xg, yg, integ, iter_f!, reinit_f!)

plot(xg,yg,basin',seriestype=:heatmap)
