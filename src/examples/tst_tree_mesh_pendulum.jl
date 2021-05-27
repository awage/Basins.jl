using Revise
using DifferentialEquations
using DynamicalSystems
using Basins
using Plots
using Printf
using RegionTrees

# Equations of motion:
function forced_pendulum(u, p, t)
    @inbounds begin
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] - sin(u[1])+ F*cos(omega*t)
    return SVector{2}(du1, du2)
    end
end

# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    if integrator.u[1] < 0
        integrator.u[1] += 2*π
    else
        integrator.u[1] -= 2*π
    end
end

condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)

cb = DiscreteCallback(condition,affect!)

#d, F ,w
F = 1.66
ω = 1.
d=0.2
p=[d, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(forced_pendulum,rand(2),(0.0,20.0), p)
integ  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-8, save_everystep=false, callback=cb)

xres=120
yres=120

# range for forced pend
xg = range(-pi,pi,length=xres)
yg = range(-2.,4.,length=yres)

@time bsn=Basins.basins_map2D_tree(xg, yg, integ; T=2π/ω, r_init=0.2, r_max=0.1)

function evaluate(basin,v)
    return findleaf(basin,v).data
end


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


plt = plot(xlim=(-pi, pi), ylim=(-2., 4.), legend=nothing)
x = range(-pi, pi, length=1000)
y = range(-2., 4., length=1000)
heatmap!(plt, x, y, (x, y) -> evaluate(bsn.basin, SVector(x, y)), fill=true)
