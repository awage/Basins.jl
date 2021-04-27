using Revise
using Plots
using Basins
using DynamicalSystems

# struct MagPendulum{T<:AbstractFloat}
#     magnets::Vector{SVector{2, T}}
# end
#
# function (m::MagPendulum)(u, p, t)
#     x, y, vx, vy = u
#     γ, d, α, ω = p
#     dx, dy = vx, vy
#     dvx, dvy = @. -ω^2*(x, y) - α*(vx, vy)
#     for ma in m.magnets
#         δx, δy = (x - ma[1]), (y - ma[2])
#         D = sqrt(δx^2 + δy^2 + d^2)
#         dvx -= γ*(x - ma[1])/D^5
#         dvy -= γ*(y - ma[2])/D^5
#     end
#     return SVector(dx, dy, dvx, dvy)
# end
#
# function mag_pendulum(u = [sincos(0.12553*2π)..., 0, 0];
#     γ = 1.0, d = 0.3, α = 0.2, ω = 0.5, N = 3)
#     m = MagPendulum([SVector(cos(2π*i/N), sin(2π*i/N)) for i in 1:N])
#     p = [γ, d, α, ω]
#     ds = ContinuousDynamicalSystem(m, u, p)
# end


# https://cdn.ima.org.uk/wp/wp-content/uploads/2020/03/Chaos-in-the-Magnetic-Pendulum-from-MT-April-2020.pdf
#ds = mag_pendulum(γ=1, d=0.5, α=0.175, ω=1., N=4)

ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=0.5, N=4)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-9)
xg=range(-2,2,length=500)
yg=range(-2,2,length=500)
@time bsn = Basins.basin_general_ds(xg, yg, integ; dt=1., idxs=1:2)


# Uncertainty exponent for these parameter and grid
bd = box_counting_dim(xg, yg, bsn)
α = 2 - bd

plot(xg,yg,bsn.basin',seriestype=:heatmap)
