using Revise
using DifferentialEquations


# Equations of motion:
function duffing!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
end

F = 0.201010
ω = 0.50202
#d, F ,w

p=[0.15, F, ω]
#p=[0.15, 0.2, 0.1]
df = ODEProblem(duffing!,rand(2),(0.0,1000.0), p)
integ_df  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false)

xres=50
yres=50
xg = range(-2.2,2.2,length=xres)
yg = range(-2.,2.,length=yres)

# compute basin
@time basin = draw_basin(xg, yg, integ_df; T=2*pi/ω)

# Basin entropy
@show Sb,Sbb = basin_entropy(basin, 20, 20)


# Basin entropy
@show Sb,Sbb = basin_entropy(basin, 20, 20)


# before computing wada merge we remove the attractors from the basin:
ind = findall(iseven.(basin) .== true)
[ basin[k]=bsn.basin[k]+1 for k in ind ]

@time max_dist,min_dist = wada_merge.compute_merge_dist(basin,xg,yg)
epsilon = xg[2]-xg[1]
@show v_dist = max_dist/epsilon


# Uncertainty dimension for these parameter and grid
@time D, ε, f_ε = uncertainty_dimension(xg, yg, integ_df; T=2*pi/ω, max_res=4, num_step=5)
