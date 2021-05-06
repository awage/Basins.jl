using Revise
using DifferentialEquations
using DynamicalSystems
using Basins
using Plots
using Printf

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
    if integrator.u[1] < 0
        integrator.u[1] += 2*π
    else
        integrator.u[1] -= 2*π
    end
end

condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)

cb = DiscreteCallback(cbasins_map2Daffect!)

#d, F ,w
F = 1.66
ω = 1.
d=0.2
p=[d, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(forced_pendulum,rand(2),(0.0,20.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false, callback=cb)

xres=200
yres=200

# range for forced pend
xg = range(-pi,pi,length=xres)
yg = range(-2.,4.,length=yres)

# compute basin
@time bsn = Basins.basins_map2D(xg, yg, integ_df; T=2*pi/ω)

# Basin entropy
@show Sb,Sbb = basin_entropy(bsn; eps_x=20, eps_y=20)

# Wada merge Haussdorff distances
@time max_dist,min_dist = detect_wada_merge_method(xg, yg, bsn)
epsilon = xg[2]-xg[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon

# Wada grid
W = detect_wada_grid_method(integ_df, bsn; max_iter=8)
@show W[:,end]

# Uncertainty exponent for these parameter and grid
bd = box_counting_dim(xg, yg, bsn)
α = 2 - bd

D = uncertainty_exponent(bsn, integ_df)
@show 2-D


println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")

@printf("Basin entropy %.2f \n", Sb)
@printf("Boundary Basin Entropy: %.2f\n", Sbb)
@printf("Uncertainty exponent: α= %.2f\n", α )
@printf("Box counting dim: bd= %.2f\n", bd)
@printf("Uncertainty dim estimator: d = %.2f\n", 2-D[1])
@printf("Number of basins: %d\n", bsn.Na)
@printf("Merge Method: Max fattening parameter: %.2f\n", dmax)
@printf("Wada Grid Method: W_Na = %.2f\n ", W[end,end] )

plot(xg, yg, bsn.basin', seriestype=:heatmap)
