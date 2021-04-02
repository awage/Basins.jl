
function box_counting_dim(xg, yg, basin, size=[])

    bnd = get_boundary_filt(basin)
    I1 = findall(bnd .== 1);
    v = hcat([[xg[ind[1]]; yg[ind[2]]] for ind in I1 ]...)
    v = Dataset(v')

    return generalized_dim(v; q=0)

end


function uncertainty_dimension_sample(xg, yg, basin)

    nx=length(xg)
    ny=length(yg)
    xi = xg[1]; xf = xg[end]; yi = yg[1];  yf = yg[end]
    grid_res_x=length(xg)
    grid_res_y=length(yg)

    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution

    min_ε = 5;
    max_ε = grid_res_x/20
    num_step=10
    r,c= size(basin)
    vals = unique(basin)
    S=Int16(length(vals))
    if S < 2
        return 1,0,0
    end

    r_ε = range(min_ε,max_ε,length=num_step)

    for (k,s_box) in enumerate(r_ε)
        Nb=0; Nu=0;
        completed = 0;
        s_box = Int(round(s_box))
        # Find uncertain boxes
        while completed == 0
            # Random box
            x = rand(1:(r-s_box+1))
            y = rand(1:(c-s_box+1))
            x_coor=x:x+s_box-1
            y_coor=y:y+s_box-1
            box_values = [basin[k,m] for k in x_coor, m in y_coor]
            Nb=Nb+1
            tmp_Nu = Nu + (length(unique(box_values))>1)
            if abs((Nu - tmp_Nu)/Nu) < 0.001 && Nu > r*c/s_box^2 && (length(unique(box_values))>1)
                completed = 1
            end
            Nu = tmp_Nu
        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,(xg[2]-xg[1])*s_box)
    end
    # uncertain function
    f_ε = N_u./N

    # get exponent
    @. model(x, p) = p[1]*x+p[2]
    fit = curve_fit(model, vec(log.(ε)), vec(log.(f_ε)), [2., 2.])
    D = coef(fit)

    return D[1], ε, f_ε

end
