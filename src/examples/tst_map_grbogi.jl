
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

θg=range(0,2π,length=256*2)
xg=range(-0.5,0.5,length=256)

@time bsn=Basins.basin_map(θg, xg, integ)
@time bns2=ChaosTools.basin_map(θg, xg, integ)

#plot(θg,xg,bsn.basin',seriestype=:heatmap)
# Estimation using box counting
bd = box_counting_dim(θg,xg, bsn.basin)
α = 2 - bd


# Estimation using the method of the paper.
# The results do not match exactly. But the paper talks about
# approximate results
N=[]
eps_r =  10 .^range(-7,0,length=10)

for eps in eps_r
    @show eps
    N2=0
    N1=0

    for k in 1:8192*2

        k1 = Int.(rand(1:(length(θg)/2)))
        k2 = rand(1:length(xg))

        c1 = bsn.basin[k1,k2]
        c2 = Basins.get_color_point!(bsn, integ, [θg[k1],xg[k2]+eps])
        c3 = Basins.get_color_point!(bsn, integ, [θg[k1],xg[k2]-eps])

        if length(unique([c1,Int(c2),Int(c3)]))>1
            N1 += 1
        end
        N2 += 1
    end
    @show N1/N2
    push!(N,N1/N2)
end

using LsqFit
# get exponent
@. model(x, p) = p[1]*x+p[2]
fit = curve_fit(model, vec(log10.(eps_r)), vec(log10.(N)), [2., 2.])
@show D = coef(fit)
@show α
plot( vec(log10.(eps_r)), vec(log10.(N)))
