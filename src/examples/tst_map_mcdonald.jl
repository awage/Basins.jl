using Revise
using DynamicalSystems
using Basins
using Plots
using BenchmarkTools


# FRACTAL BASIN BOUNDARIES
# Steven W. MCDONALD, Celso GREBOGY Edward OTT and James A. YORKE
# Physica D 17 1985
function chaotic_map(dz,z, p, n)

    xn = z[1]; yn = z[2]
    λ₁=p[1];λ₂=p[2];
    dz[1]=mod(λ₁*xn,1)
    dz[2]=λ₂*yn+cos(2*π*xn)
    return
end

# dummy function to keep the initializator happy
function chaotic_map_J(J,z0, p, n)
   return
end
ds = DiscreteDynamicalSystem(chaotic_map,[1., 0.], [1.1, 1.] , chaotic_map_J)
integ  = integrator(ds)

# This map has two escapes to ±∞
function escape_function(y)
    if y > 10000
        return 1
    elseif y < -10000
        return 2
    end
    return 0
end

function  get_color(integ,x,y)
    reinit!(integ,[x,y])
    while escape_function(integ.u[2]) == 0
        step!(integ)
        integ.t > 2000 && break
    end
    return escape_function(integ.u[2]);
end


xg=range(0.,1.,length=600)
yg=range(-2.,2.,length=600)
integ.p[1] = 3
integ.p[2] = 1.5
basin_t = zeros(length(xg),length(yg))

# Compute escape basin
for (k,x) in enumerate(xg), (n,y) in enumerate(yg)
    basin_t[k,n]=get_color(integ,x,y);
end

@show bd = box_counting_dim(xg, yg, basin_t)


# Estimation using the method of the paper.
# The results do not match exactly. But the paper talks about
# approximate results
N=[]
eps_r =  10 .^range(-8,-1,length=10)

for eps in eps_r
    @show eps
    N2=0
    N1=0

    for k in 1:8192*20

        k1 = rand(1:length(xg))
        k2 = rand(1:length(yg))

        c1 = get_color(integ, xg[k1],yg[k2])
        c2 = get_color(integ, xg[k1]+eps,yg[k2])
        c3 = get_color(integ, xg[k1]-eps, yg[k2])

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
D = coef(fit)
@show 2-D[1]

plot(xg,yg,basin_t',seriestype=:heatmap)
