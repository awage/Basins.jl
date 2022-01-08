using Revise
using Plots
using DynamicalSystems
using OrdinaryDiffEq

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
xg = range(-3.,3.,length=500)
yg = range(-3.,3.,length=500)
df_de = (;reltol = 1e-10, abstol = 1e-10, alg = Vern9())
pmap = poincaremap(ds, (3, 0.); rootkw = (xrtol = 1e-8, atol = 1e-8), diffeq = df_de)
bsn, att = basins_of_attraction((xg, yg), pmap; diffeq = df_de, mx_chk_loc_att = 100)
plot(xg, yg, bsn',seriestype = :heatmap)
