using Revise
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
# F=0.2771812080536913; ω=0.1;  # smooth boundary
#ω=0.1617;F = 0.395
#ω=0.3;F = 0.1
ds = ContinuousDynamicalSystem(duffing, rand(2), [0.15, F, ω])
xg = yg = range(-2.2,2.2,length=400)
grid = (xg,yg)
@time bsn, att = basins_of_attraction(grid, ds; T=2*pi/ω, diffeq = (;reltol = 1e-9, alg = Vern9()))

# Basin entropy
@show Sb,Sbb = basin_entropy(bsn)

# Basins fractal test
test_res, Sbb = basins_fractal_test(bsn; ε = 5, Ntotal = 10000)

# Wada merge Haussdorff distances
@show max_dist,min_dist = detect_wada_merge_method(bsn)

# Wada grid
W = detect_wada_grid_method(grid, ds; attractors = att, basins = bsn, max_iter=8, T = 2*pi/ω, diffeq = default_diffeq)
@show W[:,end]

f,e,α = uncertainty_exponent(bsn)
@show α
