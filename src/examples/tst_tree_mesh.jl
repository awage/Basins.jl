using Revise
using Plots
using Basins
using DynamicalSystems
using RegionTrees

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


ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [3] , newton_map_J)
integ  = integrator(ds)
res=10
xg=range(-1.,1.,length=res)
yg=range(-1.,1.,length=res)

@time bsn=Basins.basins_map2D_tree(xg, yg, integ)

plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), legend=nothing)
col=Plots.palette(:tab20)
for leaf in allleaves(bsn.basin)
    v = hcat(collect(vertices(leaf.boundary))...)
    #plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]])
    #@show leaf.data
    plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]],fill=(1,col[leaf.data]), linecolor=:black, lw=0.5)

end
plt
