using Revise
using Basins
using DifferentialEquations
using Plots

# Equations of motion: E. Ott, et al. I Physica D 76 (1994) 384-410
function forced_particle!(du, u, p, t)
γ=0.05  ; x̄ = 1.9  ; f₀=2.3  ; ω =3.5
x₀=1. ; y₀=0.;

    x=u[1]; dx = u[2]; y=u[3]; dy = u[4];

    du[1] = dx
    du[2] = -γ*dx -(-4*x*(1-x^2) + y^2) +  f₀*sin(ω*t)*x₀
    du[3] = dy
    du[4] = -γ*dy -(2*y*(x+x̄)) +  f₀*sin(ω*t)*y₀
end


df = ODEProblem(forced_particle!,rand(4),(0.0,20.0), [])
integ  = init(df, alg=Tsit5(); reltol=1e-9, save_everystep=false)

xres=300
yres=300

# range for forced pend
xg = range(0.,1.2,length=xres)
yg = range(0.,1.2,length=yres)

ω =3.5

# compute basin
#@time bsn = Basins.basin_general_ds(xg, yg, integ; dt=2*pi/ω, idxs=[1,3])

#plot(xg,yg,bsn.basin', seriestype=:heatmap)

function escape_function(u)

    if abs(u[3]) > 20 && abs(u[4])>100 && u[3]*u[4]>0
        @show integ.t
        return 2
    elseif abs(u[3]) < 1e-8 && abs(u[4])<1e-9 && u[3]*u[4]<0
        return 1
    end
    return 0
end

basin_t = zeros(length(xg),length(yg))
for (k,x) in enumerate(xg), (n,y) in enumerate(yg)
    reinit!(integ,[x,0,y,0])
    step!(integ, 2π/ω*10,true)
    while escape_function(integ.u) == 0
        step!(integ, 2π/ω,true)
        integ.t > 200 && break
    end

    basin_t[k,n]=escape_function(integ.u);
end

using JLD

@save "riddle.jld" basin_t
