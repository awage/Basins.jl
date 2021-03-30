Basins.jl
=========

This Julia package computes basins of attraction of dynamical systems in the phase plane and also
several metrics over the basins. The algorithm computes the basin without prior knowledge of the attractors.
This package depends heavily on the package DynamicalSystems and DifferentialEquations for the definitions of
the dynamical systems. However it is possible to define custom integrators.

The package provides the following metrics:

- Uncertainty exponent
- Basin entropy
- Wada detection
- Basin stability


## Computing the basins of attraction

The technique used to compute the basin of attraction is described in ref. [1]. It
consists in tracking the trajectory on the plane and coloring the points of according to
the attractor it leads to. This technique is very efficient for 2D basins.

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
xg=range(-4,4,length=200)
yg=range(-4,4,length=200)
basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
```

The keyword arguments are:
* `T` : the period of the stroboscopic map.
* `idxs` : the indices of the variable to track on the plane. By default the initial conditions of other variables are set to zero.


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
p= [F, G, a,b]
ds = ContinuousDynamicalSystem(lorenz84, rand(3), p)
integ  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
```

Once the integrator has been set, the Poincaré map can defined on a plane:

```jl
xg=range(-1.,1.,length=200)
yg=range(-1.5,1.5,length=200)
basin = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2)
```

The keyword arguments are:
* `plane` : A `Tuple{Int, <: Number}`, like `(j, r)` : the plane is defined
  as when the `j` variable of the system equals the value `r`. It can also be
  a vector of length `D+1`. The first `D` elements of the
  vector correspond to ``\\mathbf{a}`` while the last element is ``b``.
* `idxs`: the indices of the variable to track on the plane. By default the initial conditions of other variables are set to zero.


### Custom differential equations and low level functions.

Supose we want to define a custom ODE and compute the basin of attraction on a defined
Poincaré map:

```jl
using DifferentialEquations
using DynamicalSystems
using Basins


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

d=0.15; F=0.2; ω = 0.5848
p= [d, F, ω]
ds = ContinuousDynamicalSystem(duffing, rand(2), p)
integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)

xres=200
yres=200

xg = range(-2.2,2.2,length=xres)
yg = range(-2.2,2.2,length=yres)

iter_f! = (integ) -> step!(integ, 2π/ω, true)
reinit_f! =  (integ,y) ->  reinit!(integ, [y...])

basin = draw_basin(xg, yg, integ, iter_f!, reinit_f!)
```

The following anonymous functions are important:
* iter_f! : defines a function that iterates the system one step on the map.
* reinit_f! : sets the initial conditions on the map. Remember that only the
initial conditions on the map must be set.


## Compute the Basin Entropy

The Basin Entropy is a measure of the impredcitibility of the system when considering the initial conditions, see ref. [2].

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
basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

eps_x = 20; eps_y = 20;  
Sb,Sbb = basin_entropy(basin, eps_x, eps_y)
```
The arguments of `basin_entropy` are:
* `basin` : The basin computed on a grid.
* `eps_x`, `eps_y` : size of the window that samples the basin to compute the entropy.


## Compute the uncertainty exponent of a basin of attraction

The [uncertainty exponent](https://en.wikipedia.org/wiki/Uncertainty_exponent) and is conected to the [box-counting dimension](https://en.wikipedia.org/wiki/Box-counting_dimension). For a given resolution of the original basin, a sampling of the basin is done until the the fraction of uncertain boxes converges. The process is repeated for different box sizes and then the exponent is estimated.


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

bd = box_counting_dim(xg, yg, basin)
# uncertainty exponent
ue = 2-bd
```




## Detect the property of Wada

The [Wada property](https://en.wikipedia.org/wiki/Lakes_of_Wada) in basins of attraction is an amazing feature of some basins.
It is not trivial at all to demonstrate rigurously this property. There are however computational approaches that gives hints
about the presence of this property in a basin of attraction. One of the fastest approach is the [Merging Method](https://doi.org/10.1038/s41598-018-28119-0). The algorithm gives the maximum and minimum Haussdorff distances between merged basins. A good rule of thumb to discard the
Wada property is to check if the maximum distance is large in comparison to the resolution of the basin, i.e., if the number of pixel
is large.

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
basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

max_dist,min_dist = wada_merge_dist(basin,xg,yg)
# grid resolution
epsilon = xg[2]-xg[1]
# if dmax is large then this is not Wada
@show dmax = max_dist/epsilon
@show dmin = min_dist/epsilon
```
