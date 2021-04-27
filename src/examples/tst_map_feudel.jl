using Revise
using DynamicalSystems
using Basins
using Plots



# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
function chaotic_map(dz, z, p, n)
    xn = z[1]
    θ = z[2]
    a = p[1]
    ω = (sqrt(5.0) - 1.0) / 2.0
    r = p[2]
    f(x) = r * x * (1.0 - x)
    Mn(n) = reduce(∘, fill(f, n))
    M = Mn(3)
    dz[1] = M(xn) + a * cos(2 * π * θ)
    dz[2] = mod(θ + ω, 1.0)
    return
end

# dummy function to keep the initializator happy
function chaotic_map_J(J, z0, p, n)
    return
end
ds = DiscreteDynamicalSystem(
    chaotic_map,
    [1.0, 0.0],
    [0.0015, 3.833],
    chaotic_map_J,
)
integ = integrator(ds)


θ = range(0.0, 1.0, length = 250)
xg = range(0.0, 1.0, length = 250)
integ.p[2] = 3.833
integ.p[1] = 0.0015

@time bsn = Basins.basin_map(xg, θ, integ)
#@time bsn=Basins.basin_general_ds(xg, θ, integ; dt=1, Ncheck=10)

@show bd = box_counting_dim(xg, θ, bsn)
#@show α = 2 - bd

D = uncertainty_exponent(bsn, integ)
@show 2-D

plot(θ, xg, bsn.basin, seriestype = :heatmap)

#
# u0=[0.5, 0.];
# v = trajectory(ds,1000,u0)
#
# plot!(v[:,2],v[:,1],seriestype=:scatter)
#
# att1=Array{Float64,1}[]
# att2=Array{Float64,1}[]
# att3=Array{Float64,1}[]
#
# for v in bsn.attractors
# if v[1] == 1.
#     push!(att1,[v[2],v[3]])
# elseif v[1] == 2.
#     push!(att2,[v[2],v[3]])
# #elseif v[1] == 3.
# #    push!(att3,[v[2],v[3]])
# end
#
# end
#
# # att1 = Dataset(att1);
# # att2 = Dataset(att2);
# # att3 = Dataset(att3);
#
#
#
#
# function escape_function(u,u1,u2,u3)
#
#     m1 = minimum([norm(v - u) for v in u1])
#     m2 = minimum([norm(v - u) for v in u2])
#     #m3 = minimum([norm(v - u) for v in u3])
#     m3=1
#     if m1 < 5e-3
#         return 1
#     elseif m2 < 5e-3
#         return 2
#     elseif m3 < 5e-3
#         return 3
#     end
#     return 0
# end
#
#
# basin_t = zeros(length(xg),length(θ))
# for (k,x) in enumerate(xg), (n,y) in enumerate(θ)
#     reinit!(integ,[x,y])
#
#     while escape_function(integ.u,att1,att2,0) == 0
#         step!(integ)
#         integ.t > 200 && break
#     end
#
#     @show basin_t[k,n]=escape_function(integ.u,att1,att2,att3);
# end
#
# plot(θ,xg, basin_t, seriestype=:heatmap)
#


#
# # Estimation using box counting
# ind  = findall(iseven.(bsn.basin) .== true)
# basin_test = deepcopy(bsn.basin)
# [basin_test[k] =basin_test[k]+1 for k in ind ]
# bd = box_counting_dim(xg,θ, basin_test)
# @show α = 2 - bd
# Value in the paper is α=0.07
# But the basins are not reproductible. At least not with this method
#plot(θ,xg,bsn.basin, seriestype=:heatmap)
#
#
# r=3.7
# f(x)=r*x*(1-x)
# Mn(n) = reduce(∘, fill(f, n))
# fnn = Mn(4)
# x=0:0.01:1
# plot(x,fnn.(x))
# plot!([0,1],[0,1])
#
