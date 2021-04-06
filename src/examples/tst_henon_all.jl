using Revise
using DifferentialEquations
using DynamicalSystems
using Basins
using Plots
using Printf


ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
integ_df  = integrator(ds)

xres=200
yres=200

# range for forced pend
xg = range(-2.,2.,length=xres)
yg = range(-2.,2.,length=yres)

# compute basin
@time basin = basin_discrete_map(xg, yg, integ_df)

# Basin entropy
@show Sb,Sbb = basin_entropy(basin, 20, 20)

# Wada merge Haussdorff distances
@time max_dist,min_dist = wada_merge_dist(basin,xg,yg)
epsilon = xg[2]-xg[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon

# Wada grid
W = compute_wada_W(xg, yg, integ_df, basin; T=1, max_iter=8)
W=W./sum(W[:,1])
@show W[:,end]

# Uncertainty exponent for these parameter and grid
@time D, ε, f_ε = uncertainty_dimension_sample(xg, yg, basin)

plot(xg,yg,basin', seriestype=:heatmap)

println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")

@printf("Basin entropy %.2f \n", Sb)
@printf("Boundary Basin Entropy: %.2f\n", Sbb)
@printf("Uncertainty exponent: α= %.2f\n", D )
@printf("Number of basins: %d\n", length(unique(basin)))
@printf("Merge Method: Max fattening parameter: %.2f\n", dmax)
@printf("Wada Grid Method: W_Na = %.2f\n ", W[end,end] )


plot(xg,yg,basin', seriestype=:heatmap)
