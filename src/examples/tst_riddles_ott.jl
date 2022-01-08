using Revise
using OrdinaryDiffEq
using DynamicalSystems
using Plots

ds = System.riddled_basins(zeros(4);  γ=0.05, x̄ = 1.9, f₀=2.3, ω =3.5, x₀=1, y₀=0)
integ  = init(ds; reltol=1e-9, save_everystep=false)
res=300
ω =3.5
xg = range(-2,2,length=res)
yg = range(0.,2,length=res)
default_diffeq = (reltol = 1e-9,  alg = Vern9())
basins, att = basins_of_attraction((xg, yg), df; T = 2π/ω, diffeq = default_diffeq, horizon_limit = 10)

#
# # compute basin
# #@time bsn = Basins.basin_general_ds(xg, yg, integ; dt=2*pi/ω, idxs=[1,3])
#
# #plot(xg,yg,bsn.basin', seriestype=:heatmap)
#
# function escape_function(u)
#
#     if abs(u[3]) > 20 && abs(u[4])>100 && u[3]*u[4]>0
#         @show integ.t
#         return 2
#     elseif abs(u[3]) < 1e-8 && abs(u[4])<1e-9 && u[3]*u[4]<0
#         return 1
#     end
#     return 0
# end
#
# basin_t = zeros(length(xg),length(yg))
# for (k,x) in enumerate(xg), (n,y) in enumerate(yg)
#     reinit!(integ,[x,0,y,0])
#     step!(integ, 2π/ω*10,true)
#     while escape_function(integ.u) == 0
#         step!(integ, 2π/ω,true)
#         integ.t > 200 && break
#     end
#
#     basin_t[k,n]=escape_function(integ.u);
# end
#
# using JLD
#
# @save "riddle.jld" basin_t
