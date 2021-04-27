
using Revise
using Plots
using Basins
using DynamicalSystems


# C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
function grebogi_map(dz,z, p, n)
    θ = z[1]; x = z[2]
    J₀=0.3; a=1.32; b=0.9;
    dz[1]= θ + a*sin(2*θ) - b*sin(4*θ) -x*sin(θ)
    dz[1] = mod(dz[1],2π) # to avoid problems with attracto at θ=π
    dz[2]=-J₀*cos(θ)
    return
end

# dummy function to keep the initializator happy
function grebogi_map_J(J,z0, p, n)
   return
end
ds = DiscreteDynamicalSystem(grebogi_map,[1., -1.], [] , grebogi_map_J)
integ  = integrator(ds)

θg=range(0,2π,length=250)
xg=range(-0.5,0.5,length=250)

@time bsn=Basins.basin_map(θg, xg, integ)
#@time bns2=ChaosTools.basin_map(θg, xg, integ)

#plot(θg,xg,bsn.basin',seriestype=:heatmap)

# Estimation using box counting
@show bd = box_counting_dim(θg,xg, bsn.basin)
α = 2 - bd

# Estimation using original method
D = uncertainty_exponent(bsn, integ)
@show 2-D
