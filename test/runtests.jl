using Basins
using Test
using DynamicalSystems

@testset "Test basin_stroboscopic_map" begin
    ω=0.5
    ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
    integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
    xg=range(-2,2,length=100)
    yg=range(-2,2,length=100)
    bsn_nfo=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
    @test length(unique(bsn_nfo.basin)) == 3
    @test count(bsn_nfo.basin .== 1) == 3285
    @test count(bsn_nfo.basin .== 2) == 3285
    @test count(bsn_nfo.basin .== 3) == 3430
end

@testset "Test basin_poincare_map" begin
    b=0.1665
    ds = Systems.thomas_cyclical(b = b)
    integ=integrator(ds, reltol=1e-9)
    xg=range(-6.,6.,length=100)
    yg=range(-6.,6.,length=100)
    bsn_nfo = basin_poincare_map(xg, yg, integ; plane=(3, 0.), idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8))
    @test length(unique(bsn_nfo.basin)) == 3
    @test count(bsn_nfo.basin .== 1) == 4604
    @test count(bsn_nfo.basin .== 2) == 2813
    @test count(bsn_nfo.basin .== 3) == 2583
end

@testset "Test basin_discrete_map" begin
    ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
    integ_df  = integrator(ds)
    xg = range(-2.,2.,length=100)
    yg = range(-2.,2.,length=100)
    bsn_nfo = basin_discrete_map(xg, yg, integ_df)
    @test length(unique(bsn_nfo.basin)) == 2
    @test count(bsn_nfo.basin .== 1) == 4269
    @test count(bsn_nfo.basin .== -1) == 5731
end


@testset "Test basin_entropy" begin
    ω=0.5
    ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
    integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
    xg=range(-2,2,length=100)
    yg=range(-2,2,length=100)
    bsn_nfo = basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
    Sb,Sbb = basin_entropy(bsn_nfo.basin, 20, 20)
    @test (trunc(Sb;digits=3) == 0.853)
    @test (trunc(Sbb;digits=3) == 0.889)
end

@testset "Test wada_merge_dist" begin
    ω=0.5
    ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
    integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
    xg=range(-2,2,length=100)
    yg=range(-2,2,length=100)
    bsn_nfo=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
    # Wada merge Haussdorff distances
    max_dist,min_dist = wada_merge_dist(bsn_nfo.basin,xg,yg)
    epsilon = xg[2]-xg[1]
    dmax = max_dist/epsilon
    dmin = min_dist/epsilon
    @test (round(dmax) == 11.)
    @test (round(dmin) == 11.)
end


@testset "Test basin_stability" begin
    ω=0.5
    ds = Systems.magnetic_pendulum(γ=1, d=0.3, α=0.2, ω=ω, N=3)
    integ = integrator(ds, u0=[0,0,0,0], reltol=1e-14)
    xg=range(-2,2,length=100)
    yg=range(-2,2,length=100)
    bsn_nfo=basin_stroboscopic_map(xg, yg, integ; T=2π/ω, idxs=1:2)
    bs = basin_stability(bsn_nfo.basin)
    @test (trunc(bs[1];digits=3) == 0.328)
    @test (trunc(bs[2];digits=3) == 0.328)
    @test (trunc(bs[3];digits=3) == 0.343)
end
