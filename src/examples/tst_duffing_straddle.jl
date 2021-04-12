using Revise
using Basins
using Plots
using DynamicalSystems
using DifferentialEquations

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end


# Riddled basin, problematic
ω=0.1617
F = 0.395
ds = ContinuousDynamicalSystem(duffing, rand(2), [0.15, F, ω])
integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
xg = range(-2.2,2.2,length=100)
yg = range(-2.2,2.2,length=100)


@time bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)



sa,sb = compute_saddle(bsn, integ_df, [1], [2],100)





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
@time bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω)

Na = length(unique(bsn.basin))

sa,sb = compute_saddle(bsn, integ_df, [1], [2,3],1000)

plot(xg,yg,bsn.basin', seriestype=:heatmap)
s = Dataset(sa)
plot!(s[:,1],s[:,2],seriestype=:scatter)
