using Revise
using Plots
using DifferentialEquations
using Basins


# Equations of motion:
function forced_pendulum!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] - sin(u[1])+ F*cos(omega*t)
end


function affect!(integrator)
    if integrator.u[1] < 0
        integrator.u[1] += 2*pi
    else
        integrator.u[1] -= 2*pi
    end
end

condition(u,t,integrator) = (integrator.u[1] < -pi  || integrator.u[1] > pi)
cb = DiscreteCallback(condition,affect!)

#d, F ,w
d=0.2; F = 1.66; ω = 1.;
p=[d, F, ω]
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false, callback=cb)

# range for forced pend
xg = range(-pi,pi,length=200)
yg = range(-2.,4.,length=200)
@time basin = basin_stroboscopic_map(xg,yg,integ_df; T=2*pi/ω, idxs=1:2)
plot(xg,yg,basin', seriestype=:heatmap)
