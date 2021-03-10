using Plots
using DifferentialEquations
using Revise

include("compute_ba.jl")
import .ba_routine
#
# include("compute_wada_basins.jl")
# import .wada_merge



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
    #  integrator.u[1] = rem2pi(integrator.u[1] , RoundNearest) # wrap in -pi, pi
  #println("EVENT ", integrator.u[1])
end

condition(u,t,integrator) = (integrator.u[1] < -pi  || integrator.u[1] > pi)

cb = DiscreteCallback(condition,affect!)
#cb = ContinuousCallback(condition,affect!)

#d, F ,w
F = 1.66
ω = 1.
d=0.2
p=[d, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false, callback=cb)
# step!(integ_df, 2*pi, true)
#
# sol = solve(df,Tsit5(),callback=cb)
# plot(sol)

xres=400
yres=400

# range for forced pend
xg = range(-pi,pi,length=xres)
yg = range(-2.,4.,length=yres)


@time bsn = ba_routine.draw_basin(xg,yg,integ_df; T=2*pi/ω)

plot(xg,yg,bsn.basin', seriestype=:heatmap)