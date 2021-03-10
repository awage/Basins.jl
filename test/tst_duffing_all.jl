using Revise
using DifferentialEquations
using Basins
using Printf
# Equations of motion:
function duffing!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
end

F = 0.4
ω = 2.
#d, F ,w

p=[0.15, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)
integ_df  = init(df, alg=Tsit5(); reltol=1e-8, save_everystep=false)

xres=200
yres=200
xg = range(-2.2,2.2,length=xres)
yg = range(-2.,2.,length=yres)

# compute basin
@time basin = draw_basin(xg, yg, integ_df; T=2*pi/ω)

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
xres=50
yres=50
nxg = range(-2.2,2.2,length=xres)
nyg = range(-2.,2.,length=yres)
@time D, ε, f_ε = uncertainty_dimension(nxg, nyg, integ_df; T=2*pi/ω, max_res=5, num_step=6)

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
