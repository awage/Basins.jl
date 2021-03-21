using Revise
using Plots
using DynamicalSystems
using DifferentialEquations
using Basins


for b in 0.152:0.0005:0.178
	ds = Systems.thomas_cyclical(b = b)
	iter_f!,integ = poincaremap(ds, (3, 0.), 20., direction=+1, idxs=[1,2], rootkw = (xrtol = 1e-8, atol = 1e-8),reltol=1e-9)
	reinit_f! = (integ,y) -> reinit!(integ,[y...,0.], t0=0.)

	xg=range(-6.,6.,length=400)
	yg=range(-6.,6.,length=400)

	@time basin=draw_basin(xg, yg, integ, iter_f!, reinit_f!)

	plot(xg,yg,basin',seriestype=:heatmap)
	savefig(string("th_b_",b,".pdf"))
end
