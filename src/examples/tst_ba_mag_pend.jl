using Revise
using Plots
using Basins
using DynamicalSystems

ds = Systems.magnetic_pendulum(γ=1, d=0.2, α=0.2, ω=0.8, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-9)
xg=range(-4,4,length=350)
yg=range(-4,4,length=350)
@time bsn = basin_general_ds(xg, yg, integ; dt=1., idxs=1:2)
plot(xg,yg,bsn.basin',seriestype=:heatmap)
