using Revise
using DifferentialEquations
using DynamicalSystems
using Basins
using Plots
using Printf


ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
integ_df  = integrator(ds)

xres=300
yres=300

# range for forced pend
xg = range(-2.,2.,length=xres)
yg = range(-2.,2.,length=yres)

# compute basin
@time basin = Basins.basins_map2D(xg, yg, integ_df)

# Basin entropy
@show Sb,Sbb = basin_entropy(bsn.basin; eps_x=20, eps_y=20)

# Wada merge Haussdorff distances
@time max_dist,min_dist = detect_wada_merge_method(xg, yg, bsn.basin)
epsilon = xg[2]-xg[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon

# Wada grid
W = detect_wada_grid_method(integ_df, bsn; max_iter=8)
@show W[:,end]

# Uncertainty exponent for these parameter and grid
bd = box_counting_dim(xg, yg, bsn.basin)
α = 2 - bd

println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")

@printf("Basin entropy %.2f \n", Sb)
@printf("Boundary Basin Entropy: %.2f\n", Sbb)
@printf("Uncertainty exponent: α= %.2f\n", α )
@printf("Number of basins: %d\n", length(unique(bsn.basin)))
@printf("Merge Method: Max fattening parameter: %.2f\n", dmax)
@printf("Wada Grid Method: W_Na = %.2f\n ", W[end,end] )


plot(xg,yg,bsn.basin', seriestype=:heatmap)
