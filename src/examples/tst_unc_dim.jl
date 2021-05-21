using Revise
using Plots
using Basins
using DynamicalSystems

function newton_map(dz,z, p, n)
    f(x) = x^p[1]-1
    df(x)= p[1]*x^(p[1]-1)
    z1 = z[1] + im*z[2]
    dz1 = f(z1)/df(z1)
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end

# dummy function to keep the initializator happy
function newton_map_J(J,z0, p, n)
   return
end



# C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
function grebogi_map(dz,z, p, n)
    θ = z[1]; x = z[2]
    J₀=0.3; a=1.32; b=0.9;
    dz[1]= θ + a*sin(2*θ) - b*sin(4*θ) -x*sin(θ)
    dz[1] = mod(dz[1],2π) # to avoid problems with attracto at θ=π
    dz[2]=-J₀*cos(θ)
    return
end

# dummy function to keep the initializator happy
function grebogi_map_J(J,z0, p, n)
   return
end


function test_NM()
    ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [3] , newton_map_J)
    integ  = integrator(ds)

    r_res = range(100,700,step=100)
    bd = zeros(size(r_res))
    α =  zeros(size(r_res))
    st =  zeros(size(r_res))

    for (k,res) in enumerate(r_res)

        xg=range(-1.5,1.5,length=res)
        yg=range(-1.5,1.5,length=res)

        @time bsn=Basins.basins_map2D(xg, yg, integ)

        @show α[k],αe,f,e = uncertainty_exponent(bsn, integ)
        @show bd[k] = box_counting_dim(xg,yg, bsn)
        @show st[k],f,e= Basins.static_estimate(xg,yg,bsn)

    end

    return α,bd,st
end

function tst_grb()
    ds = DiscreteDynamicalSystem(grebogi_map,[1., -1.], [] , grebogi_map_J)
    integ  = integrator(ds)

    θg=range(0,2π,length=250)
    xg=range(-0.5,0.5,length=250)

    @time bsn=Basins.basins_map2D(θg, xg, integ)
    #@time bns2=ChaosTools.basin_map(θg, xg, integ)

    #plot(θg,xg,bsn.basin',seriestype=:heatmap)

    # Estimation using box counting
    @show bd = box_counting_dim(θg,xg, bsn.basin)
    α = 2 - bd

    # Estimation using original method
    D,De,e,f = uncertainty_exponent(bsn, integ)
    @show 2-D

    D2,e,f = Basins.static_estimate(θg,xg,bsn)
    @show 2-D2

end






ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [3] , newton_map_J)
integ  = integrator(ds)

res = 200
xg=range(-1.5,1.5,length=res)
yg=range(-1.5,1.5,length=res)

@time bsn=Basins.basins_map2D(xg, yg, integ)

@show α,αe,e,f = uncertainty_exponent(bsn, integ)
@show bd = box_counting_dim(xg,yg, bsn)
#@show st,f,e= Basins.static_estimate(xg,yg,bsn)
