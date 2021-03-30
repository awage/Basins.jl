using Basins
using Test
using DynamicalSystems

@testset "Basins.jl" begin

    println("Test basin_stroboscopic_map")
    ω=0.5
    ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
    integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
    xg=range(-2,2,length=100)
    yg=range(-2,2,length=100)
    basin=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)

    @test length(unique(basin)) == 3
    @test count(basin .== 1) == 3285
    @test count(basin .== 2) == 3285
    @test count(basin .== 3) == 3430

end
