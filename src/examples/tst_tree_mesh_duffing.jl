using Revise
using Plots
using Basins
using DynamicalSystems
using RegionTrees


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

F=0.3818791946308725; ω= 0.1966442953020134
#F=0.2771812080536913; ω=0.1;
#ω=0.1617;F = 0.395
#ω=0.3;F = 0.1
ds = ContinuousDynamicalSystem(duffing, rand(2), [0.15, F, ω])
integ_df  = integrator(ds; reltol=1e-8, save_everystep=false)
xg = range(-2.2,2.2,length=300)
yg = range(-2.2,2.2,length=300)

@time bsn=Basins.basins_map2D_tree(xg, yg, integ; r_init=0.1, r_max=0.001)

function evaluate(basin,v)
    return findleaf(basin,v).data
end

plt = plot(xlim=(-2.2, 2.2), ylim=(-2.2, 2.2), legend=nothing)
x = range(-2.2, 2.2, length=1000)
y = range(-2.2, 2.2, length=1000)
heatmap!(plt, x, y, (x, y) -> evaluate(bsn.basin, SVector(x, y)), fill=true)


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


# plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), legend=nothing)
# col=Plots.palette(:tab20)
# for leaf in allleaves(bsn.basin)
#     v = hcat(collect(vertices(leaf.boundary))...)
#     #plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]])
#     #@show leaf.data
#     plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]],fill=(1,col[mod(leaf.data,20)]), linecolor=:black, lw=0)
#
# end
# plt
