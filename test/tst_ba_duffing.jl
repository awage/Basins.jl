using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end



F=0.2
ω = 0.5848
p= [0.15, F, ω]
ds = ContinuousDynamicalSystem(duffing, rand(2), p)
integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)

xres=200
yres=200

xg = range(-2.2,2.2,length=xres)
yg = range(-2.2,2.2,length=yres)

@time basin = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)
plot(xg,yg,basin', seriestype=:heatmap)
