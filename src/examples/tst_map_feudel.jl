using Revise
using DynamicalSystems
using Basins
using Plots




# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
function chaotic_map(dz,z, p, n)

    xn = z[1]; θ = z[2]
    a=p[1]; ω = (sqrt(5.)-1.)/2.; r=p[2]
    f(x)=r*x*(1. - x)
    M(x)=f(f(f(x)))
    dz[1]=M(xn)+a*cos(2*π*θ)
    dz[2]=mod(θ + ω, 1.)
    return
end

# dummy function to keep the initializator happy
function chaotic_map_J(J,z0, p, n)
   return
end
ds = DiscreteDynamicalSystem(chaotic_map,[1., 0.], [0.0015, 3.833] , chaotic_map_J)
integ  = integrator(ds)


θ=range(0.,1.,length=200)
xg=range(0.,1.,length=200)
integ.p[2]=3.846
integ.p[1] = 0.0024

@time bsn=Basins.basin_map(xg, θ, integ)
#@time bsn=Basins.basin_general_ds(xg, θ, integ; dt=1, Ncheck=10)

plot(θ,xg,bsn.basin, seriestype=:heatmap)


u0=[0.5, 0.];
v = trajectory(ds,1000,u0)

plot!(v[:,2],v[:,1],seriestype=:scatter)

att1=Array{Float64,1}[]
att2=Array{Float64,1}[]
att3=Array{Float64,1}[]

for v in bsn.attractors
if v[1] == 1.
    push!(att1,[v[2],v[3]])
elseif v[1] == 2.
    push!(att2,[v[2],v[3]])
#elseif v[1] == 3.
#    push!(att3,[v[2],v[3]])
end

end

# att1 = Dataset(att1);
# att2 = Dataset(att2);
# att3 = Dataset(att3);




function escape_function(u,u1,u2,u3)

    m1 = minimum([norm(v - u) for v in u1])
    m2 = minimum([norm(v - u) for v in u2])
    #m3 = minimum([norm(v - u) for v in u3])
    m3=1
    if m1 < 5e-3
        return 1
    elseif m2 < 5e-3
        return 2
    elseif m3 < 5e-3
        return 3
    end
    return 0
end


basin_t = zeros(length(xg),length(θ))
for (k,x) in enumerate(xg), (n,y) in enumerate(θ)
    reinit!(integ,[x,y])

    while escape_function(integ.u,att1,att2,0) == 0
        step!(integ)
        integ.t > 200 && break
    end

    @show basin_t[k,n]=escape_function(integ.u,att1,att2,att3);
end

plot(θ,xg, basin_t, seriestype=:heatmap)



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


#sa,sb = compute_saddle(integ, bsn, [1], [2,3]; N=1000)

# r=3.846
# f(x)=r*x*(1-x)
# M(x)=f(f(f(x)))
# n=0:0.01:1
# plot(n,M.(n)); plot!([0,1],[0,1])




# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
# function chaotic_map_2(z)
#     xn = z[1]; θ = z[2]
#     dz=[0. 0.]
#     a=0.0024; ω = (sqrt(5)-1)/2; r=3.846
#     f(x)=r*x*(1-x)
#     M(x)=f(f(f(x)))
#
#     dz[1]=M(xn)+a*cos(2*π*θ)
#     dz[2]=mod(θ + ω,1)
#     return dz
# end


#
# u0=[0.46,0.45]
# N=1000
# v=zeros(N,2)
# u=u0
# for k=1:N
#
#     v[k,:]=chaotic_map_2(u)
#     u=v[k,:]
# end
#
# v2=trajectory(ds, 1000, u0)
# reinit!(integ,u0)
# v3 =zeros(N,2)
# for k=1:N
#     step!(integ)
#     v3[k,:]=integ.u
# end
#
#
# plot(v[:,2],v[:,1],seriestype=:scatter)
#
# plot!(v2[:,2],v2[:,1],seriestype=:scatter)
# plot!(v3[:,2],v3[:,1],seriestype=:scatter)
