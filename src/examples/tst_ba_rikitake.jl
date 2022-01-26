using Revise
using Plots
using DynamicalSystems
using OrdinaryDiffEq

μ=0.47
ds = Systems.rikitake(μ = μ, α = 1.0)
xg = yg = range(-2.5,2.5,length = 100)
df_de = (;reltol = 1e-10, abstol = 1e-10, alg = Vern9(), maxiters = Int(1e12))
pmap = poincaremap(ds, (3, 0.); rootkw = (xrtol = 1e-8, atol = 1e-8), diffeq = df_de)
bsn, att = basins_of_attraction((xg, yg), ds; mx_chk_loc_att = 100)
plot(xg, yg, bsn',seriestype = :heatmap)

# This system is problematic: very very long transient (t > 2000 sometimes) that can be mistaken with attractors.


 delete!(att, 3) # remove "wrong" attractor
 delete!(att, 4) # remove "wrong" attractor

 xg = yg = range(-2.5,2.5,length = 400)
 bsn, att = basins_of_attraction((xg, yg), ds; attractors = att)
 plot(xg, yg, bsn',seriestype = :heatmap)
