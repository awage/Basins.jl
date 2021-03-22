# Basins

#[![Build Status](https://github.com/awage/Basins.jl/workflows/CI/badge.svg)](https://github.com/awage/Basins.jl/actions)


Basins.jl
=========

This Julia package computes basins of attraction of dynamical systems in the phase plane and also
several metrics over the basins. The algorithm computes the basin without prior knowledge of the attractors.

The package provides the following metrics:

- Uncertainty exponent
- Basin entropy
- Wada detection
- Basin stability


## Computing the basins of attraction


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



* iter_f! : defines a function that iterates the system one step on the map.
* reinit_f! : sets the initial conditions on the map. Remember that only the
initial conditions on the map must be set.
