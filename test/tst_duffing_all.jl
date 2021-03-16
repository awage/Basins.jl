using Revise
using Basins
using Printf
using Plots
using DynamicalSystems
using DifferentialEquations

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

F=0.2
ω = 0.5848
p= [0.15, F, ω]
ds = ContinuousDynamicalSystem(duffing, rand(2), p)
integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)


xres=200
yres=200
xg = range(-2.2,2.2,length=xres)
yg = range(-2.,2.,length=yres)

# compute basin
iter_f! = (x) -> step!(x, 2*pi/ω, true)
reinit_f! = (x, y) -> reinit!(x, y, t0=0, erase_sol=true,reinit_callbacks=true)
@time basin=draw_basin(xg, yg, integ_df, iter_f!, reinit_f!)


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
