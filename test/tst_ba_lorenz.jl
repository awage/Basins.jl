using Revise
using Plots
using Basins
using DifferentialEquations
using Printf
using DynamicalSystems
using Roots


function poincaremap(integ, planecrossing, Tmax, j, rootkw)
    f = (t) -> planecrossing(integ(t))
	ti = integ.t
    #data = _initialize_output(integ.u, j)
    #Ttr != 0 && step!(integ, Ttr)

    # Check if initial condition is already on the plane
    side = planecrossing(integ.u)
    if side == 0
        #push!(data, integ.u[j])
		dat = integ.u[j]
        step!(integ)
        side = planecrossing(integ.u)
		return dat
    end

    #while integ.t < tfinal + Ttr
        while side < 0
            (integ.t - ti) > Tmax && break
            step!(integ)
            side = planecrossing(integ.u)
			#@show integ.u, integ.t

        end
        while side ≥ 0
            (integ.t - ti) > Tmax && break
            step!(integ)
            side = planecrossing(integ.u)
			#@show integ.u, integ.t
        end

		# Did not found the crossing. Tmax reached.
		# Maybe a fixed point outside the plane
		 if (integ.t - ti) > Tmax
		 	#@show integ.u,integ.t
		 	return integ.u[j]
		 end

        # I am now guaranteed to have `t` in negative and `tprev` in positive
        tcross = Roots.find_zero(f, (integ.tprev, integ.t), Roots.A42(); rootkw...)
        ucross = integ(tcross)
        #push!(data, ucross[j])
    #end
    return ucross[j]
end


# Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model
# Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
@inline @inbounds function lorenz84(u, p, t)
    F = p[1]; G = p[2]; a = p[3]; b = p[4];
	x = u[1]; y = u[2]; z = u[3];
    dx = -y^2 -z^2 -a*x + a*F
    dy = x*y - y - b*x*z +G
	dz = b*x*y + x*z -z
    return SVector{3}(dx, dy, dz)
end


F=6.846
G=1.287
a=0.25
b=4.
p= [F, G, a,b]
ds = ContinuousDynamicalSystem(lorenz84, rand(3), p)
integ  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)

#integ = integrator(ds, u0=[0,0,0])
plane = (3, 0.)
direction = -1
planecrossing = PlaneCrossing(plane, direction > 0)
rootkw = (xrtol = 1e-8, atol = 1e-8)
idxs=[1,2]
iter_f! = (integ) -> poincaremap(integ, planecrossing, 20,  SVector{length(idxs), Int}(idxs...), rootkw)
#iter_f! = (integ) -> step!(integ, 0.01)
reinit_f! = (integ,y) -> reinit!(integ,[y...,0.], t0=0)

xg=range(-1.,1.,length=400)
yg=range(-1.5,1.6,length=400)
#zg=range(0,50,length=100)

reinit_f!(integ,[0.,0.])
#u = Dataset([ iter_f!(integ) for k in 1:2000])
#plot(u[:, 1], u[:, 2], seriestype=:scatter)


@time basin=draw_basin(xg, yg, integ, iter_f!, reinit_f!)

plot(xg,yg,basin',seriestype=:heatmap)
