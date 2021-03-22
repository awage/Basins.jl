using Revise
using Plots
using Basins
using DifferentialEquations
using DynamicalSystems



function newton_map(dz,z, p, n)
    f(x) = x^p[1]-1
    df(x)= p[1]*x^(p[1]-1)
    z1 = z[1] + im*z[2]
    dz1 = f(z1)/df(z1)
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end

# dummy function to keep the initializator happy
function newton_map_J(J,z0, p, n)
   return
end

ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [10] , newton_map_J)
integ  = integrator(ds)

xg=range(-1.5,1.5,length=200)
yg=range(-1.5,1.5,length=200)

@time basin=basin_discrete_map(xg, yg, integ)
plot(xg,yg,basin',seriestype=:heatmap)
