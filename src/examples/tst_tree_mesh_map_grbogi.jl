using Revise
using Plots
using Basins
using DynamicalSystems
using RegionTrees

# C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
ds = Systems.grebogi_map(rand(2))
integ  = integrator(ds)

θg=range(0,2π,length=1000)
xg=range(-1,1,length=1000)

@time bsn=Basins.basins_map2D_tree(θg, xg, integ; r_init=0.2, r_max=0.01)

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
    ε = reverse(wd * 2 .^ range(0,k_mx-1,step=1))
    return ε,N[1:k_mx]
end


ε,N = compute_frac_dim(bsn)

elog = log10.(1 ./ ε)
Nlog = log10.(N)
D = linear_region(elog[4:end],Nlog[4:end])
@show D


plt = plot(xlim=(θg[1], θg[end]), ylim=(xg[1], xg[end]), legend=nothing)
heatmap!(plt, θg, xg, (x, y) -> evaluate(bsn.basin, SVector(x, y)), fill=true)

#
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
