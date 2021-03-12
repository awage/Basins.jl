using Plots
using Basins
using DifferentialEquations
using Printf


ODE_SOLVER = Vern7
#ODE_SOLVER = Tsit5
#HAM_SOLVER = DPRKN6
HAM_SOLVER = Yoshida6

function henon_heiles_de(dv,v,p,t)
    x = v[1];
    px = v[2];
    y = v[3];
    py = v[4];
    ∂H∂x = x+2*x*y;
    ∂H∂y = y-y^2+x^2;
    ∂H∂px = px;
    ∂H∂py = py;

	dv[1] = ∂H∂px;
	dv[2] = -∂H∂x;
	dv[3] = ∂H∂py;
	dv[4] = -∂H∂y;
end


function salida(sol)
    a=sol[1,end]
    b=sol[3,end]
    if b>1
        sal=1;
    elseif b<0 && a<-1
        sal=2;
    elseif b<0 && a>1
        sal=3;
    else sal=0;
    end
    return sal
end




function get_exit(Ei, y0, py, x0)
    Esqrt = 2. *Ei-y0^2+2/3*y0^3-py^2;
    if Esqrt < 0
        return 0
    end
    px=sqrt(Esqrt);
    vi = [x0,px,y0,py];
    tspan=(0.0,400)
    # Define callback to halt solver
    condition(u,t,integrator) = (u[1]^2+u[3]^2)>100
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prob = ODEProblem(henon_heiles_de,vi,tspan,callback=cb)
    #sol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12)
    sol = solve(prob,ODE_SOLVER(),reltol=1e-8,abstol=1e-8)
    return salida(sol)
end


Ei=0.25
x0=0.
ry0=range(-0.65,1.1,length=300)
rpy=range(-0.75,0.75,length=300)

@time basin = [  get_exit(Ei, y0, py, x0) for y0=ry0, py=rpy]

plot(ry0,rpy,basin', seriestype=:heatmap)

# Basin entropy
@show Sb,Sbb = basin_entropy(basin, 20, 20)

# Wada merge Haussdorff distances
@time max_dist,min_dist = wada_merge_dist(basin,ry0,rpy)
epsilon = ry0[2]-ry0[1]
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon


println("---------------")
println("---------------")
println("Basin Report: ")
println("---------------")
println("---------------")

@printf("Basin entropy %.2f \n", Sb)
@printf("Boundary Basin Entropy: %.2f\n", Sbb)
@printf("Number of basins: %d\n", length(unique(basin)))
@printf("Merge Method: Max fattening parameter: %.2f\n", dmax)




using Roots


function poincaremap(integ, planecrossing, Ttr, j, rootkw)
    f = (t) -> planecrossing(integ(t))
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
            #integ.t > tfinal + Ttr && break
            step!(integ)
            side = planecrossing(integ.u)
        end
        while side ≥ 0
            #integ.t > tfinal + Ttr && break
            step!(integ)
            side = planecrossing(integ.u)
        end
        #integ.t > tfinal + Ttr && break

        # I am now guaranteed to have `t` in negative and `tprev` in positive
        tcross = Roots.find_zero(f, (integ.tprev, integ.t), Roots.A42(); rootkw...)
        ucross = integ(tcross)
        #push!(data, ucross[j])
    #end
    return ucross[j]
end



hh = Systems.henonheiles()
plane = (1, 0)
u0 = [0.0, -0.25, 0.42081, 0.0]
integ = integrator(hh, u0=u0)
direction = -1
planecrossing = PlaneCrossing(plane, direction > 0)
rootkw = (xrtol = 1e-6, atol = 1e-6)
reinit!(integ, [0.0, -0.25, 0.42081, 0.0])
idxs=[2,4]
#poincarmap(integ, planecrossing, 0, SVector{length(idxs), Int}(idxs...), rootkw)
f_iter = (integ) -> poincaremap(integ, planecrossing, 0,  SVector{length(idxs), Int}(idxs...), rootkw)

u = Dataset([ f_iter(integ) for k in 1:2000])

plot(u[:, 1], u[:, 2], seriestype=:scatter)
