using Revise
using OrdinaryDiffEq
using DynamicalSystems
using Basins
using Plots
using Printf

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
d=0.2
ds = Systems.forced_pendulum([0.,0.]; ω = ω, f = F, d = d);
# range for forced pend
res = 300
xg = range(-pi, pi,length = res)
yg = range(-2., 4.,length = res)
grid = (xg,yg)
default_diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
bsn, att = basins_of_attraction(grid, ds; T = 2*pi/ω, diffeq = default_diffeq)

# Basin entropy
@show Sb,Sbb = basin_entropy(bsn)

# Basins fractal test
test_res, Sbb = basins_fractal_test(bsn; ε = 5, Ntotal = 10000)

# Wada merge Haussdorff distances
@time max_dist,min_dist = detect_wada_merge_method(xg, yg, bsn)
epsilon = xg[2]-xg[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon

# Wada grid
W = detect_wada_grid_method(grid, ds; attractors = att, basins = bsn, max_iter=8, T = 2*pi/ω, diffeq = default_diffeq)
@show W[:,end]

# Uncertainty exponent for these parameter and grid
ε, N_ε ,α = uncertainty_exponent(bsn)
@show α


println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")
#
# @printf("Basin entropy %.2f \n", Sb)
# @printf("Boundary Basin Entropy: %.2f\n", Sbb)
# @printf("Uncertainty exponent: α= %.2f\n", α )
# @printf("Box counting dim: bd= %.2f\n", bd)
# @printf("Uncertainty dim estimator: d = %.2f\n", 2-D[1])
# @printf("Number of basins: %d\n", bsn.Na)
# @printf("Merge Method: Max fattening parameter: %.2f\n", dmax)
# @printf("Wada Grid Method: W_Na = %.2f\n ", W[end,end] )

plot(xg, yg, basins', seriestype=:heatmap)
