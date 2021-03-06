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
ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [3] , newton_map_J)
integ  = integrator(ds)

a=zeros(1,8)
b=zeros(1,8)
c=zeros(1,8)
d=zeros(1,8)

# for (k,res) in enumerate(range(100,800,step=100))
#
#     xg=range(-1.,1.,length=res)
#     yg=range(-1.,1.,length=res)
#
#     @time bsn=Basins.basins_map2D(xg, yg, integ)
#
#     a[k],De,e,f=uncertainty_exponent(bsn,integ)
#     b[k],e,f=Basins.static_estimate(xg,yg,bsn)
#     c[k],e,f=Basins.ball_estimate(xg,yg,bsn)
#     d[k]=2-box_counting_dim(xg,yg,bsn)
# end

xg=range(-1.,1.,length=300)
yg=range(-1.,1.,length=300)

@time bsn=Basins.basins_map2D(xg, yg, integ)


b=zeros(1,40)
c=zeros(1,40)

precision = 1e-5
for k in 1:40
        #a[k],De,e,f=uncertainty_exponent(bsn,integ)
        b[k],e,f=Basins.static_estimate(xg,yg,bsn; precision=precision)
        c[k],e,f=Basins.ball_estimate(xg,yg,bsn; precision=precision)
        #d[k]=2-box_counting_dim(xg,yg,bsn)
end


using Statistics
println(mean(b), "  ", mean(c))
#plot(xg,yg,bsn.basin',seriestype=:heatmap)
plot(b')
plot!(c')


#sa,sb = compute_saddle(integ, bsn, [1], [2,3]; N=100)

#s=Dataset(sa)
#plot!(s[:,1],s[:,2],seriestype=:scatter)
