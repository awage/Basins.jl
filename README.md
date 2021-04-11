Basins.jl
=========

This Julia package computes basins of attraction of dynamical systems in the phase plane and also
several metrics over the basins. The algorithm computes the basin without prior knowledge of the attractors.

This package depends heavily on the package DynamicalSystems and DifferentialEquations for the definitions of the dynamical systems. However it is possible to define custom integrators.

The package provides the following metrics:

- Uncertainty exponent
- Basin entropy
- Wada detection
- Basin stability


## 1 - Computing the basins of attraction

The technique used to compute the basin of attraction is described in ref. [1]. It consists in tracking the trajectory on the plane and coloring the points of according to the attractor it leads to. This technique is very efficient for 2D basins.

The algorithm gives back a matrix with the attractor numbered from 1 to N. If an attractor exists outside the defined grid of if the trajectory escapes, this initial condition is labelled -1. It may happens for example if there is a fixed point not on the Poincaré map.

### Usage

First define a dynamical system on the plane, for example with a stroboscopic map or Poincaré section. For example we can set up an dynamical system with a stroboscopic map defined:

```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
```

Now we define the grid of ICs that we want to analyze and launch the procedure:

```jl
xg=range(-2,2,length=200)
yg=range(-2,2,length=200)
bsn=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
```

The keyword arguments are:
* `T` : the period of the stroboscopic map.
* `idxs` : the indices of the variable to track on the plane. By default the initial conditions of other variables are set to zero.

The function returns a structure `bsn` with several fields of interests:
* `bsn.basin` is a matrix that contains the information of the basins of attraction. The attractors are numbered from 1 to N and each element
correspond to an initial condition on the grid.
* `bsn.xg` and `bsn.yg` are the grid vectors.
* `bsn.attractors` is a collection of vectors with the location of the attractors found.

Now we can plot the nice result of the computation:

```jl
using Plots
plot(xg,yg,bsn.basin', seriestype=:heatmap)

```

https://i.imgur.com/EBWw1GK.png

Another example with a Poincaré map:
```jl
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

F=6.846; G=1.287; a=0.25; b=4.;
ds = ContinuousDynamicalSystem(lorenz84, rand(3), [F, G, a,b])
integ  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
```

Once the integrator has been set, the Poincaré map can defined on a plane:

```jl
xg=range(-1.,1.,length=200)
yg=range(-1.5,1.5,length=200)
bsn = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2)
```

The keyword arguments are:
* `plane` : A `Tuple{Int, <: Number}`, like `(j, r)` : the plane is defined
  as when the `j` variable of the system equals the value `r`. It can also be
  a vector of length `D+1`. The first `D` elements of the
  vector correspond to ``\\mathbf{a}`` while the last element is ``b``.
* `idxs`: the indices of the variable to track on the plane. By default the initial conditions of other variables are set to zero.


## 2 - Custom differential equations and low level functions.

Supose we want to define a custom ODE and compute the basin of attraction on a defined
Poincaré map:

```jl
using DifferentialEquations
using Basins

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

d=0.15; F=0.2; ω = 0.5848
ds = ContinuousDynamicalSystem(duffing, rand(2), [d, F, ω])
integ = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
xg = range(-2.2,2.2,length=200)
yg = range(-2.2,2.2,length=200)

iter_f! = (integ) -> step!(integ, 2π/ω, true)
reinit_f! =  (integ,y) ->  reinit!(integ, [y...])

bsn = draw_basin(xg, yg, integ, iter_f!, reinit_f!)
```

The following anonymous functions are important:
* iter_f! : defines a function that iterates the system one step on the map.
* reinit_f! : sets the initial conditions on the map. Remember that only the
initial conditions on the map must be set.


## 3 - Compute the Basin Entropy

The [Basin Entropy](https://doi.org/10.1007/978-3-319-68109-2_2) is a measure of the impredictability of the basin of attraction of a dynamical system. An important feature of the basins of attraction is that for a value above log(2) we can say that the basin is fractalized.

### Usage

Once the basin of attraction has been computed, the computing the Basin Entropy is easy:

```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
bsn=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

Sb,Sbb = basin_entropy(bsn.basin; eps_x=20, eps_y=20)
```
The arguments of `basin_entropy` are:
* `basin` : The basin computed on a grid.
* `eps_x`, `eps_y` : size of the window that samples the basin to compute the entropy.


## 4 - Compute the uncertainty exponent of a basin of attraction

The [uncertainty exponent](https://en.wikipedia.org/wiki/Uncertainty_exponent) is conected to the [box-counting dimension](https://en.wikipedia.org/wiki/Box-counting_dimension). For a given resolution of the original basin, a sampling of the basin is done until the the fraction of uncertain boxes converges. The process is repeated for different box sizes and then the exponent is estimated.


### Usage


```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
bsn=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

bd = box_counting_dim(xg, yg, bsn.basin)

# uncertainty exponent is the dimension of the plane minus the box-couting dimension
ue = 2-bd
```




## 5 - Detection of the property of Wada

### 5.1 - Merge Method

The [Wada property](https://en.wikipedia.org/wiki/Lakes_of_Wada) in basins of attraction is an amazing feature of some basins. It is not trivial at all to demonstrate rigurously this property. There are however computational approaches that gives hints about the presence of this property in a basin of attraction. One of the fastest approach is the [Merging Method](https://doi.org/10.1038/s41598-018-28119-0). The algorithm gives the maximum and minimum Haussdorff distances between merged basins. A good rule of thumb to discard the Wada property is to check if the maximum distance is large in comparison to the resolution of the basin, i.e., if the number of pixel is large.

Notice that the algorithm gives an answer for a particular choice of the grid. It is not an accurate method.

### Usage

```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
bsn=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

max_dist,min_dist = detect_wada_merge_method(xg, yg, bsn.basin)
# grid resolution
epsilon = xg[2]-xg[1]
# if dmax is large then this is not Wada
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon
```

### 5.2 - Grid Method

Another method available and much more accurate is the [Grid Method](https://doi.org/10.1038/srep16579). It divides the grid and scrutinize the boundary to test if all the attractors are present in every point of the boundary. It may be very long to get an answer since the number of points to test duplicates at each step. The algorithm returns a vector with the proportion of boxes with 1 to N attractor. For example if the vector W[N] is above 0.95 we have all the initial boxes in the boundary on the grid with N attractors. It is therefore a strong evidence that we have a Wada boundary.  


### Usage

```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
bsn=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

@show W = detect_wada_grid_method(integ, bsn_nfo; max_iter=10)
```

The algorithm returns:
* `W` contains a vector with the proportion of boxes in the boundary of `k` attractor. A good criterion to decide if the boundary is Wada is to look at `W[N]` with N the number of attractors. If this number is above 0.95 we can conclude that the boundary is Wada.  


## 6 - Computation of the Basin Stability

The Basin Stability [6] measures the relative sizes of the basin. Larger basin are considered more stable since a small perturbation or error in the initial conditions is less likely to change the attractor.

### Usage

```jl
using DynamicalSystems
using Basins
ω=0.5
ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

@show basin_stability(bsn.basin)
```


## References

[1] H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012

[2] A. Daza, A. Wagemakers, B. Georgeot, D. Guéry-Odelin and M. A. F. Sanjuán, Basin entropy: a new tool to analyze uncertainty in dynamical systems, Sci. Rep., 6, 31416 (2016).

[3] C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983

[4] A. Daza, A. Wagemakers and M. A. F. Sanjuán, Ascertaining when a basin is Wada: the merging method, Sci. Rep., 8, 9954 (2018).

[5] A. Daza, A. Wagemakers, M. A. F. Sanjuán and J. A. Yorke, Testing for Basins of Wada, Sci. Rep., 5, 16579 (2015).

[6] P. Menck, J. Heitzig, N. Marwan et al. How basin stability complements the linear-stability paradigm. Nature Phys 9, 89–92 (2013).
