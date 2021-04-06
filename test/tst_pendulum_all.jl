using Revise
using DifferentialEquations
using Basins
using Printf

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
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false, callback=cb)

xres=200
yres=200

# range for forced pend
xg = range(-pi,pi,length=xres)
yg = range(-2.,4.,length=yres)

# compute basin
@time basin = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω)

# Basin entropy
@show Sb,Sbb = basin_entropy(basin, 20, 20)

# Wada merge Haussdorff distances
@time max_dist,min_dist = wada_merge_dist(basin,xg,yg)
epsilon = xg[2]-xg[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon

# Wada grid
W = compute_wada_W(xg, yg, integ_df, basin; T=2*pi/ω, max_iter=8)
W=W./sum(W[:,1])
@show W[:,end]

# Uncertainty exponent for these parameter and grid
D, ε, f_ε = uncertainty_dimension_sample(xg, yg, basin)

bd = box_counting_dim(xg, yg, basin)

plot(xg,yg,basin', seriestype=:heatmap)

println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")

@printf("Basin entropy %.2f \n", Sb)
@printf("Boundary Basin Entropy: %.2f\n", Sbb)
@printf("Uncertainty exponent: α= %.2f\n", D )
@printf("Box counting dim: bd= %.2f\n", bd)

@printf("Number of basins: %d\n", length(unique(basin)))
@printf("Merge Method: Max fattening parameter: %.2f\n", dmax)
@printf("Wada Grid Method: W_Na = %.2f\n ", W[end,end] )


plot(xg,yg,basin', seriestype=:heatmap)
