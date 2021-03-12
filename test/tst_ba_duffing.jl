using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins

#
# # Equations of motion:
# function duffing!(du, u, p, t)
#     d = p[1]; F = p[2]; omega = p[3]
#     du[1] = u[2]
#     du[2] = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
#
# end
#
#
#
# # Cuenca rara, a veces se detectan 3 y otras 2 atractores.
# #F=0.419192
# #ω = 0.10
#
# F=0.2
# ω = 0.5848
#
# #d, F ,w
# p=  [0.15, F, ω]
# df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)
#
# integ_df  = init(df, alg=Tsit5(); reltol=1e-8, save_everystep=false)


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

F=0.2
ω = 0.5848

p= [0.15, F, ω]

ds = ContinuousDynamicalSystem(duffing, rand(2), p)

integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
#integ_df  = init(df, alg=Tsit5(); reltol=1e-9, save_everystep=false)



xres=200
yres=200

xg = range(-2.2,2.2,length=xres)
yg = range(-2.2,2.2,length=yres)

@time basin = draw_basin(xg, yg, integ_df; T=2*pi/ω)

plot(xg,yg,basin', seriestype=:heatmap)
