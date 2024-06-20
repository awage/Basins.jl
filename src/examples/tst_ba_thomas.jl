using Revise
using Plots
using DynamicalSystems
using OrdinaryDiffEq

### First version using high-level function.
b=0.1665
#b=0.22
ds = Systems.thomas_cyclical(b = b)
xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)

pmap = poincaremap(ds, (3, 0.), 1e6;  diffeq = (;reltol = 1e-9, alg = Vern9()))
@time bsn, att = basins_of_attraction((xg, yg), pmap)

plot(xg, yg, bsn', seriestype = :heatmap)
