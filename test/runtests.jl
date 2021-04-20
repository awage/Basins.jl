using Basins
using DynamicalSystems
using Test
using DifferentialEquations

@testset "All tests" begin

@testset "Test basin_stroboscopic_map" begin
    ω=1.; F = 0.2
    ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
    integ_df  = integrator(ds; alg=AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-8, save_everystep=false)
    xg = range(-2.2,2.2,length=100)
    yg = range(-2.2,2.2,length=100)
    @time bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)
    @test length(unique(bsn.basin))/2 == 2
    @test count(bsn.basin .== 3) == 5376
    @test count(bsn.basin .== 5) == 4622

end

@testset "Test basin_poincare_map" begin
    b=0.1665
    ds = Systems.thomas_cyclical(b = b)
    xg=range(-6.,6.,length=100)
    yg=range(-6.,6.,length=100)
    pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9)
    bsn = basin_poincare_map(xg, yg, pmap)

    @test length(unique(bsn.basin))/2 == 3
    @test count(bsn.basin .== 3) == 4639
    @test count(bsn.basin .== 5) == 2680
    @test count(bsn.basin .== 7) == 2665
end

@testset "Test basin_discrete_map" begin
    ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
    integ_df  = integrator(ds)
    xg = range(-2.,2.,length=100)
    yg = range(-2.,2.,length=100)
    bsn_nfo = basin_discrete_map(xg, yg, integ_df)

    @test count(bsn_nfo.basin .== 3) == 4127
    @test count(bsn_nfo.basin .== -1) == 5730
end


@testset "Test basin_entropy" begin
    ω=1.; F = 0.2
    ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
    integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
    xg = range(-2.2,2.2,length=150)
    yg = range(-2.2,2.2,length=150)
    bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)
    Sb,Sbb = basin_entropy(bsn.basin; eps_x=20, eps_y=20)
    @test (trunc(Sb;digits=3) == 0.663)
    @test (trunc(Sbb;digits=3) == 0.663)
end

@testset "Test Wada detection" begin
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

    F = 1.66;   ω = 1.;    d=0.2
    df = ODEProblem(forced_pendulum!,rand(2),(0.0,20.0), [d, F, ω])
    integ  = init(df, alg=AutoTsit5(Rosenbrock23()); reltol=1e-9, save_everystep=false, callback=cb)
    xg = range(-pi,pi,length=100); yg = range(-2.,4.,length=100)
    bsn = basin_stroboscopic_map(xg, yg, integ; T=2*pi/ω)

    # Wada merge Haussdorff distances
    # First remove attractors
    ind  = findall(iseven.(bsn.basin) .== true)
    basin_test = deepcopy(bsn.basin)
    [basin_test[k] =basin_test[k]+1 for k in ind ]
    max_dist,min_dist = detect_wada_merge_method(xg,yg,basin_test)
    epsilon = xg[2]-xg[1]
    dmax = max_dist/epsilon
    dmin = min_dist/epsilon
    @test (round(dmax) == 6.)
    @test (round(dmin) == 2.)

    W = detect_wada_grid_method(integ, bsn; max_iter=5)
    @test (trunc(W[3];digits=2) ==  0.93)

end


@testset "Test basin_stability" begin

    ω=1.; F = 0.2
    ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
    integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
    xg = range(-2.2,2.2,length=150)
    yg = range(-2.2,2.2,length=150)
    bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)

    bs = basin_stability(bsn.basin)
    @test (trunc(bs[1];digits=3) == 0.492)
    @test (trunc(bs[2];digits=3) == 0.507)
end


@testset "Test box_counting_dimension" begin

    ω=1.; F = 0.2
    ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
    integ_df  = integrator(ds; alg=Tsit5(),  reltol=1e-8, save_everystep=false)
    xg = range(-2.2,2.2,length=150)
    yg = range(-2.2,2.2,length=150)
    bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)
    ind  = findall(iseven.(bsn.basin) .== true)
    basin_test = deepcopy(bsn.basin)
    [basin_test[k] =basin_test[k]+1 for k in ind ]
    bd = box_counting_dim(xg, yg, basin_test)
    @test (trunc(bd;digits=2) == 1.9)
end

@testset "Test compute saddle straddle" begin
    v = [-0.7024272244361562 0.7358943669783878;
    0.6884108286512226 0.2883122496731336;
    0.12248777086662276 0.17266855751576238;
    0.21416001598855067 0.1918853810113003;
    -0.13345893435577877 -0.10181462673868892;
    -0.24612685358870395 0.40609732169904256;
    0.5443525330629784 0.3637084477432412;
    -0.8543699001914774 -0.5504191054878991;
    -0.81424974452034 0.3876573128398122;
    -0.3759850591885575 0.3029466252309591;
    1.2695931242637941 0.5649763379043485;
    0.12319548967885363 0.1988038193210813;
    0.08993268394324054 0.08763876246458675;
    0.5998194677317075 0.42267589676742806;
    -0.7471048204006432 0.4173581434012174;
    -0.3371249654066878 0.41666798624746265;
    0.7345532189737 0.3144499987176077;
    0.24865664622871408 0.28838520561505615;
    -0.7047154958073086 -0.5213541656840109;
    -0.9056118879414881 0.33965597519541424;
    -0.4547496618520956 0.21744361580215943;
    -0.2921283747093624 0.4214569617780614;
    0.6291128784123545 0.3716326220739068;
    -0.6010738837992531 -0.43440268456863923;
    -0.9652688621193279 0.24069456477396253;
    -0.5501113728424021 0.15418201057091677;
    -0.4754568961525241 0.10081275951313767;
    -0.46112021536727277 0.09967906030053533;
    -0.4391298368793252 0.11974265732897998;
    -0.3662749654329448 0.20493692475101746;
    -0.4603329215220823 0.6299851440826713;
    0.41313450130864204 0.28437379172110344;
    -0.7299720705099347 -0.5197447840442747;
    -0.890354901563265 0.3447067895447532;
    -0.44687447282979653 0.22139280314316853;
    -0.2966913365698902 0.45089092769960826;
    0.5501212667337514 0.3597852144857481;
    -0.8200621499512287 -0.5412846282080438;
    -0.8349858140921245 0.3758250755769921;
    -0.39618264893265825 0.2735219498521972;
    -0.9437331293087026 0.7574461204326174;
    0.9293611309012172 -0.15957232746111116;
    1.14599339117127 0.4412985174267645;
    1.2396517289691749 0.5188989178649657;
    0.8959919160009866 0.7090718967680737;
    -1.3461532547680735 0.008447716620129747;
    -0.6398916478847173 0.6224442163197923;
    0.8276869977268546 0.06600919221071469;
    1.0607570462732214 0.6160653104393968;
    -0.20258432736829904 -0.10525350252623214;
    -0.3349774815559638 0.19748108608150314;
    -0.663004994554249 0.724099207141831;
    0.6512663688101122 0.31465908699556383;
    -0.19176348994343262 -0.1154502893224026;
    -0.33183913666489595 0.1997476599045356;
    -0.7339332858876773 0.7450732722276884;
    0.7126484478736383 0.2607963013619859;
    0.3565945556966931 0.376716105526132;
    -1.1424583292686683 -0.5181306982978956;
    -0.6691798553091624 0.48524940958270907;
    -0.7750640627426556 0.7447635013295555;
    0.7617537229889098 0.18145588387180991;
    0.7785640313451194 0.6364791644564665;
    -1.3274874780888364 -0.13525661761177976;
    -0.6135443463613546 0.580073323763926;
    0.8897837152389886 -0.06286426899098416;
    1.1440178842431123 0.4993844802023786;
    1.071146877796164 0.6852925880897688;
    -1.0130674106682127 -0.5202505944779693;
    -0.7194964747274598 0.40984327478863253;
    -0.3337755688232476 0.43117917669428807;
    0.6995015529974392 0.3448638592935367;
    -0.09956822741216353 -0.02900414846930688;
    -0.8926894003760993 0.7683209764135135;
    0.8340505900639127 0.011479529745785512;
    1.0830174110506972 0.5892351694272199;
    0.21204311031642775 0.27421769347596214;
    -0.5644430315888005 -0.4442700978076312;
    -0.9833207542715674 0.22818523879153269;
    -0.5545537845226697 0.1593390081784191;
    -0.472231790321022 0.10490889580215614;
    -0.44780564466696554 0.11274653830575736;
    -0.39238113473288266 0.17168250535822036;
    -0.2884018330542508 0.44862621970732874;
    0.5312147137503425 0.35474079716608226;
    -0.8467395756009193 -0.5497951927484416;
    -0.8188633312213228 0.38589323804591835;
    -0.37959560576794266 0.29752006205944304;
    0.4696497379499575 0.5119045223664362;
    -1.3699491532902666 0.0061234201083442685;
    -0.7483772772691258 0.6806762370674923;
    0.8713539043026773 -0.04321198642007509;
    1.1294232129600181 0.527761159044075;
    0.8977550302213696 0.7036254552586797;
    -1.3389562103162267 -0.04653014856135522;
    -0.6203974969982126 0.6005519054452689;
    0.8447656866536019 0.030093293866878423;
    1.096993725067534 0.5832862071649934;
    0.3206215383741723 0.3693331736195564;
    -1.1173678907969715 -0.5383392012987194;
    -0.6806544022816443 0.47871124803080173;
    -0.626554882949195 0.6976130759887653;
    0.6516324418974114 0.31993125496392566;
    -0.2171788716252682 -0.13785785370964124;
    -0.43298341824677933 0.08012326524258628;
    -0.42802570899573394 0.12623349965496122;
    -0.33932641866474456 0.2438970723609983;
    -1.2365333502089975 -0.2700884633323442];

    v = v';
    ω=1.; F = 0.2
    ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
    integ_df  = integrator(ds; alg=AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-8, save_everystep=false)
    xg = range(-2.2,2.2,length=200)
    yg = range(-2.2,2.2,length=200)
    @time bsn = basin_stroboscopic_map(xg, yg, integ_df; T=2*pi/ω, idxs=1:2)
    sa,sb = compute_saddle(integ_df, bsn, [1], [2]; N=100)
    s = hcat(sa...)
    hd = Basins.haussdorff_dist(s,v)
    @test (hd < 0.4)

end

end
