using Revise
using Basins
using Plots
using DynamicalSystems
using DifferentialEquations

#
#
# ω=1.
# F = 0.2
# ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
# integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
# xg = range(-2.2,2.2,length=150)
# yg = range(-2.2,2.2,length=150)
#
# @time bsn = Basins.basins_map2D(xg, yg, integ_df; T=2*pi/ω)
#
# sa,sb = compute_saddle(integ_df, bsn, [1], [2]; N=1000)

# plot(xg,yg,bsn.basin', seriestype=:heatmap)
# s = Dataset(sa)
# plot!(s[:,1],s[:,2],seriestype=:scatter)

# io = open("myfile.txt", "w");
# for v in s
#     write(io, string(v[1], " ", v[2], "; \n"));
# end
# close(io);




# Forced pendulum

# Equations of motion:
function forced_pendulum!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] - sin(u[1])+ F*cos(omega*t)
end

# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    if integrator.u[1] < 0
        integrator.u[1] += 2*π
    else
        integrator.u[1] -= 2*π
    end
end

condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)

cb = DiscreteCallback(condition,affect!)

#d, F ,w
F = 1.66
ω = 1.
d=0.2
p=[d, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, abstol=1e-9, save_everystep=false, callback=cb)

xres=200
yres=200

# range for forced pend
xg = range(-pi,pi,length=xres)
yg = range(-2.,4.,length=yres)

# compute basin
@time bsn = Basins.basins_map2D(xg, yg, integ_df; T=2*pi/ω)

Na = length(unique(bsn.basin))

sa,sb = compute_saddle(integ_df, bsn, [1], [2,3]; N=1000)

plot(xg,yg,bsn.basin', seriestype=:heatmap)
s = Dataset(sa)
plot!(s[:,1],s[:,2],seriestype=:scatter, markercolor=:blue)
