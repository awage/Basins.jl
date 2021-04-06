using Revise
using Plots
using DynamicalSystems
using Basins


# Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model
# Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
@inline @inbounds function lorenz84(u, p, t)
    F = p[1]; G = p[2]; a = p[3]; b = p[4];
	x = u[1]; y = u[2]; z = u[3];
    dx = -y^2 -z^2 -a*x + a*F
    dy = x*y - y - b*x*z +G
	dz = b*x*y + x*z -z
    return SVector{3}(dx, dy, dz)
end


F=6.846; G=1.287; a=0.25; b=4.;
p= [F, G, a, b]
ds = ContinuousDynamicalSystem(lorenz84, rand(3), p)
integ  = integrator(ds; reltol=1e-8)

xg=range(-1.,1.,length=100)
yg=range(-1.5,1.6,length=100)

@btime basin = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2);
#plot(xg,yg,basin',seriestype=:heatmap)
