"""
    box_counting_dim(xg, yg, basin; kwargs...)
This function compute the box-counting dimension of the boundary. It is related to the uncertainty dimension.

[C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983]

## Arguments
* `basin` : the matrix containing the information of the basin.
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.

## Keyword arguments
* `kwargs...` these arguments are passed to `generalized_dim` see `ChaosTools.jl` for the docs.

"""
function box_counting_dim(xg, yg, basin; kwargs...)
    # get points coordinates in the boudary
    bnd = get_boundary_filt(basin)
    I1 = findall(bnd .== 1);
    v = Dataset{2,eltype(xg)}()
    for  ind in I1; push!(v, [xg[ind[1]] , yg[ind[2]]]); end;
    return generalized_dim(v; q=0, kwargs...)
end

function box_counting_dim(xg, yg, bsn::BasinInfo; kwargs...)
    ind  = findall(iseven.(bsn.basin) .== true)
    # Not sure it is necessary to deepcopy
    basin_test = deepcopy(bsn.basin)
    # Remove attractors from the picture
    for k in ind; basin_test[k] = basin_test[k]+1; end
    return box_counting_dim(xg, yg, basin_test; kwargs...)
end




"""
    uncertainty_exponent(bsn::BasinInfo, integ; sizes=0.)
This function estimates the uncertainty exponent of the boundary. It is related to the uncertainty dimension.

[C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983]

## Arguments
* `bsn` : structure that holds the information of the basin.
* `integ` : handle to the iterator of the dynamical system.

## Keyword arguments
* `sizes` array that specifies the scales at which the uncertainty exponent should be computed.

"""
function uncertainty_exponent(bsn::BasinInfo, integ; sizes=0.)
    xg = bsn.xg; yg = bsn.yg;
    nx=length(xg)
    ny=length(yg)
    xi = xg[1]; xf = xg[end]; yi = yg[1];  yf = yg[end]
    grid_res_x=xg[2]-xg[1]
    grid_res_y=yg[2]-yg[1]

    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution


    min_ε = 1e-6;
    max_ε = grid_res_x*5;
    num_step=15
    r,c= size(bsn.basin)
    vals = sum(isodd.(unique(bsn.basin)))
    S=Int16(vals)
    if S < 2
        return 1,0,0
    end

    r_ε = range(min_ε,max_ε,length=num_step)

    for (k,eps) in enumerate(r_ε)
        Nb=0; Nu=0; μ=0; σ²=0; M₂=0;
        completed = 0;
        # Find uncertain boxes
        while completed == 0
            k1 = rand(1:length(xg))
            k2 = rand(1:length(yg))

            c1 = bsn.basin[k1,k2]
            c2 = get_color_point!(bsn, integ, [xg[k1]+eps,yg[k2]])
            c3 = get_color_point!(bsn, integ, [xg[k1]-eps,yg[k2]])

            if length(unique([c1,Int(c2),Int(c3)]))>1
                Nu = Nu + 1
            end
            Nb += 1

            # Welford's online mean and average estimation
            M₂ = wel_var(M₂, μ, Nu/Nb, Nb)
            μ = wel_mean(μ, Nu/Nb, Nb)
            σ² = M₂/Nb

            # Stopping criterion: variance of the estimator of the mean bellow 1e-3
            if Nb > 100 && σ² < 1e-3
                completed = 1
                #@show Nu,Nb
            end

        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,eps)
    end

    # uncertain function
    f_ε = N_u./N

    # get exponent
    @. model(x, p) = p[1]*x+p[2]
    fit = curve_fit(model, vec(log10.(ε)), vec(log10.(f_ε)), [2., 2.])
    D = coef(fit)

    return D[1]
end



function wel_var(M₂, μ, xₙ, n)

    μ₂ = μ + (xₙ - μ)/n
    M₂ = M₂ + (xₙ - μ)*(xₙ - μ₂)
    return M₂
end

function wel_mean(μ, xₙ, n)
    return μ + (xₙ - μ)/n
end
