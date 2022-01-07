using DynamicalSystems
using Plots

# Equations of motion:
function forced_pendulum(u, p, t)
    @inbounds begin
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] - sin(u[1])+ F*cos(omega*t)
    return SVector{2}(du1, du2)
    end
end

# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    uu = integrator.u
    if integrator.u[1] < 0
        set_state!(integrator, SVector(uu[1] + 2π, uu[2]))
        u_modified!(integrator, true)
    else
        set_state!(integrator, SVector(uu[1] - 2π, uu[2]))
        u_modified!(integrator, true)
    end
end


condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)
cb = DiscreteCallback(condition,affect!)

#d, F ,w
F = 1.66
ω = 1.
d = 0.2
p = [d, F, ω]
# p = [0.15, 0.2, 0.1]
#p = [0.2, 1.5454545454545454, 0.6515151515151515]
#p = [0.2, 1.3636363636363635, 0.5] # Parametro para cuenca riddle.
ds = ContinuousDynamicalSystem(forced_pendulum,rand(2), p)

# range for forced pend
xg = range(-pi,pi,length = 1000)
yg = range(-15.,15.,length = 1000)
grid = (xg,yg)
default_diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
bsn, att = basins_of_attraction(grid, ds; T = 2*pi/ω, diffeq = default_diffeq)

sa,sb = compute_saddle(grid, ds; bas_A = [1], bas_B = [2,3] N=100)

plot(xg,yg,bsn.basin', seriestype=:heatmap)
s = Dataset(sa)
plot!(s[:,1],s[:,2],seriestype=:scatter, markercolor=:blue)
