using Revise
using Plots
using DynamicalSystems
using OrdinaryDiffEq

# Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model
# Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
F=6.846; G=1.287; a=0.25; b=4.;
ds = Systems.lorenz84(; F,G,a,b)
xg=range(-5.,5.,length = 500)
yg=range(-5.,5.,length = 500)
zg=range(-2.,2.,length = 10)
#pmap = poincaremap(ds, (3, 0.))
@time bsn, att = basins_of_attraction((xg, yg, zg), ds)
plot(xg,yg,bsn[:,:,5]',seriestype=:heatmap)
