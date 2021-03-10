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

#d, F ,w
F = 0.201010
ω = 0.50202
p=[0.15, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-8, save_everystep=false)

xres=200
yres=200
xg = range(-2.2,2.2,length=xres)
yg = range(-2.,2.,length=yres)

@time basin = draw_basin(xg, yg, integ_df; T=2*pi/ω)

@show Sb,Sbb = basin_entropy(basin, 20, 20)

plot(xg,yg,basin', seriestype=:heatmap)
