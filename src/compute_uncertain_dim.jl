
function uncertainty_dimension(xg, yg, integ_df; T, num_step=10, max_res=4)

    nx=length(xg)
    ny=length(yg)
    xi = xg[1]; xf = xg[end]; yi = yg[1];  yf = yg[end]
    grid_res_x=length(xg)
    grid_res_y=length(yg)

    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution

    for (k,r) in enumerate(range(1,max_res,length=num_step))

        x_n = range(xi,xf, length=Int(round(grid_res_x*r)))
        y_n = range(yi,yf, length=Int(round(grid_res_y*r)))

        @time basin = draw_basin(x_n, y_n, integ_df; T)

        r,c= size(basin)
        vals = unique(basin)
        S=Int16(length(vals))
        if S < 2
            return 1,0,0
        end
        Nb=0; Nu=0;
        eps_x=3;
        eps_y=3;
        # Find uncertain boxes
        for x = 1:eps_x:(r-eps_x+1), y = 1:eps_y:(c-eps_y+1)
                x_coor=x:x+eps_x-1
                y_coor=y:y+eps_y-1
                box_values=[basin[k,m] for k in x_coor, m in y_coor]
                Nb=Nb+1
                Nu = Nu + (length(unique(box_values))>1)
        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,x_n[2]-x_n[1])
    end

    # uncertain function
    f_ε = N_u./N

    # get exponent
    @. model(x, p) = p[1]*x+p[2]
    fit = curve_fit(model, vec(log.(ε)), vec(log.(f_ε)), [2., 2.])
    D = coef(fit)

    return D[1], ε, f_ε

end
