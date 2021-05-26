using Revise
using Basins
using Printf
using Plots
using DynamicalSystems

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

F=0.3818791946308725; ω= 0.1966442953020134
#F=0.2771812080536913; ω=0.1;
#ω=0.1617;F = 0.395
#ω=0.3;F = 0.1
ds = ContinuousDynamicalSystem(duffing, rand(2), [0.15, F, ω])
integ_df  = integrator(ds; reltol=1e-8, save_everystep=false)
xg = range(-2.2,2.2,length=300)
yg = range(-2.2,2.2,length=300)


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
