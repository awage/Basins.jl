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

@time bsn=Basins.basins_map2D_tree(xg, yg, integ; r_init=0.2, r_max=0.001)

function evaluate(basin,v)
    return findleaf(basin,v).data
end

function compute_frac_dim(bsn)

    root_n(x,n) = foldl((y,_) -> parent(y), 1:n; init=x)

    N = zeros(Int, 20)
    wd=1e3
    k_mx = 0
    for lf in allcells(bsn.basin)

        k = 1
        while !isnothing(root_n(lf,k))
            k += 1
        end
        N[k] += 1
        wd = min(wd,lf.boundary.widths[1])
        k_mx = max(k_mx,k)
    end
    @show N
    eps = reverse(wd * 2 .^ range(0,k_mx-1,step=1))
    return eps,N[1:k_mx]
end


eps,N = compute_frac_dim(bsn)

elog = log10.(1 ./ eps)
Nlog = log10.(N)
D = linear_region(elog[4:end],Nlog[4:end])
@show D




plt = plot(xlim=(xg[1], xg[end]), ylim=(yg[1], yg[end]), legend=nothing)

x = range(xg[1], xg[end], length=1000)
y = range(yg[1], yg[end], length=1000)
heatmap!(plt, x, y, (x, y) -> evaluate(bsn.basin, SVector(x, y)), fill=true)


#
# # plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), legend=nothing)
#  for leaf in allleaves(bsn.basin)
#      v = hcat(collect(vertices(leaf.boundary))...)
#      #plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]])
#      #@show leaf.data
#      plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]], linecolor=:black, lw=0.4)
#
#  end
#  plt
