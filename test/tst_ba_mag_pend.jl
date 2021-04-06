using Revise
using Plots
using Basins
using DynamicalSystems


ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-2,2,length=100)
yg=range(-2,2,length=100)
@time bsn = basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
