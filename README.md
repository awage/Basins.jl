Basins.jl
=========

This Julia package is now deprecated in favor of [ChaosTools.jl](https://github.com/JuliaDynamics/ChaosTools.jl). However some experimental stuff remains in the package: Wada detection methods and Saddles computation.


The package provides the following metrics:

1. Wada detection
2. Basin saddles
3. Examples



## 1 - Detection of the property of Wada

### 1.1 - Merge Method

The [Wada property](https://en.wikipedia.org/wiki/Lakes_of_Wada) in basins of attraction is an amazing feature of some basins. It is not trivial at all to demonstrate rigurously this property. There are however computational approaches that gives hints about the presence of this property in a basin of attraction. One of the fastest approach is the [Merging Method](https://doi.org/10.1038/s41598-018-28119-0). The algorithm gives the maximum and minimum Haussdorff distances between merged basins. A good rule of thumb to discard the Wada property is to check if the maximum distance is large in comparison to the resolution of the basin, i.e., if the number of pixel is large.

Notice that the algorithm gives an answer for a particular choice of the grid. It is not an accurate method.

### Usage

```jl
using Basins, DynamicalSystems, DifferentialEquations

# Equations of motion:
function forced_pendulum!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] - sin(u[1])+ F*cos(omega*t)
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
F = 1.66; ω = 1.; d=0.2
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), [d, F, ω])
integ = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, abstol=1e-9, save_everystep=false, callback=cb)
bsn = basins_map2D(range(-pi,pi,length=100), range(-2.,4.,length=100), integ; T=2*pi/ω)

max_dist,min_dist = detect_wada_merge_method(xg, yg, bsn)
# grid resolution
epsilon = xg[2]-xg[1]
# if dmax is large then this is not Wada
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon
```

### 1.2 - Grid Method

Another method available and much more accurate is the [Grid Method](https://doi.org/10.1038/srep16579). It divides the grid and scrutinize the boundary to test if all the attractors are present in every point of the boundary. It may be very long to get an answer since the number of points to test duplicates at each step. The algorithm returns a vector with the proportion of boxes with 1 to N attractor. For example if the vector W[N] is above 0.95 we have all the initial boxes in the boundary on the grid with N attractors. It is therefore a strong evidence that we have a Wada boundary.  


### Usage

```jl
using Basins, DynamicalSystems, DifferentialEquations

# Equations of motion:
function forced_pendulum!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] - sin(u[1])+ F*cos(omega*t)
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
F = 1.66; ω = 1.; d=0.2
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), [d, F, ω])
integ = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, abstol=1e-9, save_everystep=false, callback=cb)
bsn = basins_map2D(range(-pi,pi,length=100), range(-2.,4.,length=100), integ; T=2*pi/ω)

@show W = detect_wada_grid_method(integ, bsn; max_iter=10)
```

The algorithm returns:
* `W` contains a vector with the proportion of boxes in the boundary of `k` attractor. A good criterion to decide if the boundary is Wada is to look at `W[N]` with N the number of attractors. If this number is above 0.95 we can conclude that the boundary is Wada.  


## 2 - Computation of the saddle embedded in the boundary [EXPERIMENTAL FEATURE]

There is an invariant subset of the boundary which is invariant under the forward iteration of the dynamical system. This set is called the chaotic set, chaotic saddle or simply saddle set. It is possible to compute an approximation arbitrarily close to the saddle with the saddle straddle method. For a detailed description of the method see [1]. This method requires two `generalized basins` such that the algorithm focus on the boundary between these two sets. We divide the basins in two class such that  `bas_A ∪ bas_B = [1:N]` and `bas_A ∩ bas_B = ∅` with `N` the number of attractors.

### Usage


```jl
using Basins, DynamicalSystems, DifferentialEquations

# Equations of motion:
function forced_pendulum!(du, u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du[1] = u[2]
    du[2] = -d*u[2] - sin(u[1])+ F*cos(omega*t)
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
F = 1.66; ω = 1.; d=0.2
df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), [d, F, ω])
integ = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, abstol=1e-9, save_everystep=false, callback=cb)
bsn = basins_map2D(range(-pi,pi,length=200), range(-2.,4.,length=200), integ; T=2*pi/ω)

# sa is the left set and sb is the right set.
sa,sb = compute_saddle(integ, bsn, [1], [2,3], 1000)
s = Dataset(sa) # convert to a dataset for ploting
plot(xg,yg,bsn.basin', seriestype=:heatmap)
plot!(s[:,1],s[:,2],seriestype=:scatter, markercolor=:blue)
```


The arguments of `compute_saddle` are:
* `integ` : the matrix containing the information of the basin.
* `bsn_nfo` : structure that holds the information of the basin as well as the map function. This structure is set when the basin is first computed with `basins_map2D` or `basin_poincare_map`.
* `bas_A` : vector with the indices of the attractors that will represent the generalized basin A
* `bas_B` : vector with the indices of the attractors that will represent the generalized basin B. Notice that `bas_A ∪ bas_B = [1:N]` and `bas_A ∩ bas_B = ∅`

Keyword arguments are:
* `N` : number of points of the saddle to compute


![image](https://i.imgur.com/pQLDO0Ol.png)


## References

[1] H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 1997

[2] A. Daza, A. Wagemakers, B. Georgeot, D. Guéry-Odelin and M. A. F. Sanjuán, Basin entropy: a new tool to analyze uncertainty in dynamical systems, Sci. Rep., 6, 31416 (2016).

[3] C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983

[4] A. Daza, A. Wagemakers and M. A. F. Sanjuán, Ascertaining when a basin is Wada: the merging method, Sci. Rep., 8, 9954 (2018).

[5] A. Daza, A. Wagemakers, M. A. F. Sanjuán and J. A. Yorke, Testing for Basins of Wada, Sci. Rep., 5, 16579 (2015).
