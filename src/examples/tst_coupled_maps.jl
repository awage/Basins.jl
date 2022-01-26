using Revise
using Plots
using Basins
using DynamicalSystems


# SHRIMALI, M. D., PRASAD, A., RAMASWAMY, R., & FEUDEL, U. (2008). THE NATURE OF ATTRACTOR BASINS IN MULTISTABLE SYSTEMS. International Journal of Bifurcation and Chaos, 18(06), 1675–1688. doi:10.1142/s0218127408021269
function coupled_logistic_maps(dz,z, p, n)
    θ = z[3]; xn = z[1]; yn = z[2]
    α=3.5; β=0.1; β₂=β*(4/α-1); β₁=0.01; ω = (√5-1)/2
    dz[1]=α*(1 +β₂*cos(2π*θ))*xn*(1-xn)+β₁*(yn-xn)
    dz[2]=α*(1 +β₂*cos(2π*θ))*yn*(1-yn)+β₁*(xn-yn)
    dz[3]=mod(θ+ ω,1)
    return
end

# dummy function to keep the initializator happy
function coupled_logistic_maps_J(J,z0, p, n)
   return
end
ds = DiscreteDynamicalSystem(coupled_logistic_maps,[1., -1., 0.], [] , coupled_logistic_maps_J)
integ  = integrator(ds)

xg=range(0.,0.999999,length=100)
yg=range(0.,0.999999,length=100)

@time bsn, att = basins_of_attraction((xg, yg), ds)
plot(xg,yg,bsn',seriestype=:heatmap)

# # Estimation using box counting
# ind  = findall(iseven.(bsn.basin) .== true)
# basin_test = deepcopy(bsn.basin)
# [basin_test[k] =basin_test[k]+1 for k in ind ]
# bd = box_counting_dim(xg,yg, basin_test)
# α = 2 - bd
