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

θg=range(0,2π,length=600)
xg=range(-0.5,0.5,length=600)

@time bsn=Basins.basins_map2D_tree(θg, xg, integ; r_init=0.1, r_max=0.01)

function evaluate(basin,v)
    return findleaf(basin,v).data
end

plt = plot(xlim=(0, 2π), ylim=(-0.5, 0.5), legend=nothing)
x = range(0, 2π, length=1000)
y = range(-0.5, 0.5, length=1000)
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
