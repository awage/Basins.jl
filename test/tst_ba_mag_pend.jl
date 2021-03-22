using Revise
using Plots
using Basins
using DifferentialEquations
using DynamicalSystems


ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=300)
yg=range(-4,4,length=300)
@time basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
plot(xg,yg,basin',seriestype=:heatmap)
