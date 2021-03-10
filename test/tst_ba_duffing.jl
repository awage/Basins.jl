using Revise
using Plots
using DifferentialEquations
using Basins


# Equations of motion:
function duffing!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)

end



# Cuenca rara, a veces se detectan 3 y otras 2 atractores.
#F=0.419192
#ω = 0.10

F=0.2
ω = 0.5848

#d, F ,w
p=  [0.15, F, ω]
df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)

integ_df  = init(df, alg=Tsit5(); reltol=1e-8, save_everystep=false)

xres=200
yres=200

xg = range(-2.2,2.2,length=xres)
yg = range(-2.2,2.2,length=yres)

@time basin = draw_basin(xg, yg, integ_df; T=2*pi/ω)

plot(xg,yg,basin', seriestype=:heatmap)
