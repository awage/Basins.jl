using Revise
using Plots
using Basins
using DifferentialEquations
using DynamicalSystems

# Blowout bifurcation with non-normal parameters in population dynamicsBernard Cazelle
function franke_yabuku(dz,z, p, n)
    r=2.825; s=20.25; c1=1.2; c2=0.1;
    x = z[1]; y = z[2];
    dx = x*exp(r-s*(x+y))
    dy = c1*y/(c2+x+y)
    dz = [dx,dy]
    return
end

# dummy function to keep the initializator happy
function franke_yabuku_J(J,z0, p, n)
   return
end

ds = DiscreteDynamicalSystem(franke_yabuku,[0.1, 0.2], [] , franke_yabuku_J)
integ  = integrator(ds)

xg=range(8.,16.,length=100)
yg=range(0.5,1.,length=100)

# compute basin

@time basin=basins_map2D(xg, yg, integ)
plot(xg,yg,basin',seriestype=:heatmap)
