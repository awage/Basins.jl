
function uncertainty_dimension(xg, yg, integ_df, T)

    nx=length(xg)
    ny=length(yg)
    xi = xg[1]; xf = xg[end]; yi = yg[1];  yf = yg[end]


    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution

    for (k,r) in enumerate(range(1,4,length=10))

        x_n = range(xi,xf, length=Int(round(r*100)))
        y_n = range(yi,yf, length=Int(round(r*100)))

        @time bsn = draw_basin(x_n, y_n, integ_df; T)
        # before computing wada merge we remove the attractors from the basin (even numbers):
        ind = findall(iseven.(bsn.basin) .== true)
        [ bsn.basin[k[1],k[2]]=bsn.basin[k[1],k[2]]+1 for k in ind ]

        r,c= size(bsn.basin)
        vals = unique(bsn.basin)
        S=Int16(length(vals))
        Nb=0; Nu=0;
        eps_x=3;
        eps_y=3;
        # Find uncertain boxes
        for x = 1:eps_x:(r-eps_x+1), y = 1:eps_y:(c-eps_y+1)
                x_coor=x:x+eps_x-1
                y_coor=y:y+eps_y-1
                box_values=[bsn.basin[k,m] for k in x_coor, m in y_coor]
                Nb=Nb+1
                Nu = Nu + (length(unique(box_values))>1)
        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,x_n[2]-x_n[1])
    end

    return  ε, N_u./N


end
