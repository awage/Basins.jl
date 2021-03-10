using Plots
using DifferentialEquations
using BenchmarkTools
using Revise
using Profile

include("compute_ba.jl")
import .ba_routine

include("compute_wada_basins.jl")
import .wada_merge

# Equations of motion:
function duffing!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
end

#d, F ,w
F = 0.201010
ω = 0.50202
p=[0.15, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false)

N=5
v_dist = zeros(1,N)
#k=3
for k=1:N
    xres=k*100
    yres=k*100
    xg = range(-2.2,2.2,length=xres)
    yg = range(-2.,2.,length=yres)


    @time bsn = ba_routine.draw_basin(xg, yg, integ_df; T=2*pi/ω)
    #plot(xg,yg,bsn.basin', seriestype=:heatmap)
    # anim = @animate for ω = 0.1:0.01:0.50202
    #     df = ODEProblem(duffing!,rand(2),(0.0,1000.0), [0.15, F, ω])
    #     @time bsn=ba_routine.draw_basin(xg,yg,df; T=2*pi/ω)
    #     plot(xg,yg,bsn.basin', seriestype=:heatmap)
    # end
    #
    # gif(anim, "anim_fps15.gif", fps = 2)


    # before computing wada merge we remove the attractors from the basin:
    ind = findall(iseven.(bsn.basin) .== true)
    [ bsn.basin[k[1],k[2]]=bsn.basin[k[1],k[2]]+1 for k in ind ]
    bsn.basin = (bsn.basin .- 1)./2 # convert odd number to sequence 1:Na for the merge method

    @time max_dist,min_dist = wada_merge.compute_haussdorff_dist_n(bsn.basin,xg,yg)
    epsilon = xg[2]-xg[1]
    #@show v_dist[k]=(max_dist-min_dist)/min_dist
    @show v_dist[k]=max_dist/epsilon
end
