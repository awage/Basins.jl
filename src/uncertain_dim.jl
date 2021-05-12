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
    uncertainty_exponent(bsn::BasinInfo, integ; precision=1e-3) -> α,e,ε,f_ε
This function estimates the uncertainty exponent of the boundary. It is related to the uncertainty dimension.

[C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983]

## Arguments
* `bsn` : structure that holds the information of the basin.
* `integ` : handle to the iterator of the dynamical system.

## Keyword arguments
* `precision` variance of the estimator of the uncertainty function.

"""
function uncertainty_exponent(bsn::BasinInfo, integ; precision=1e-3)
    xg = bsn.xg; yg = bsn.yg;
    nx=length(xg)
    ny=length(yg)
    xi = xg[1]; xf = xg[end];
    grid_res_x=xg[2]-xg[1]


    # First we compute a crude estimate to have an idea at which resolution
    D,e,f = static_estimate(xg,yg,bsn);
    println("First estimate: α=",D)

    # estimate the minimum resolution to estimate a f_ε ∼ 10^-3
    @show min_ε = -3/D[1]

    min_ε = max(min_ε,-7) # limit resolution to 10^-7 to avoid numerical inconsistencies.

    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution

    # resolution in log scale
    #min_ε = -6;
    max_ε = log10(grid_res_x*5);
    num_step=10

    r_ε = 10 .^ range(min_ε,max_ε,length=num_step)

    for (k,eps) in enumerate(r_ε)
        Nb=0; Nu=0; μ=0; σ²=0; M₂=0;
        completed = 0;
        # Find uncertain boxes
        while completed == 0
            k1 = rand(1:nx)
            k2 = rand(1:ny)

            c1 = bsn.basin[k1,k2]
            c2 = get_color_point!(bsn, integ, [xg[k1]+eps,yg[k2]])
            c3 = get_color_point!(bsn, integ, [xg[k1]-eps,yg[k2]])

            if length(unique([c1,Int(c2),Int(c3)]))>1
                Nu = Nu + 1
            end
            Nb += 1

            # Welford's online average estimation and variance of the estimator
            M₂ = wel_var(M₂, μ, Nu/Nb, Nb)
            μ = wel_mean(μ, Nu/Nb, Nb)
            σ² = M₂/Nb
            # If the process is binomial, the variance of the estimator is Nu/Nb·(1-Nu/Nb)/Nb (p·q/N)
            # simulations matches
            #push!(sig,Nu/Nb*(1-Nu/Nb)/Nb)

            # Stopping criterion: variance of the estimator of the mean bellow 1e-3
            if Nu > 50 && σ² < precision
                completed = 1
                @show Nu,Nb,σ²
            end

            if Nu < 10 && Nb>10000
                # skip this value, we don't find the boundary
                # corresponds to roughly f_ε ∼ 10^-4
                @warn "Resolution may be to small for this basin, skip value"
                completed = 1
                Nu = 0
            end

        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,eps)
    end

    # uncertain function
    f_ε = N_u./N

    # remove zeros in case there are:
    ind = f_ε .> 0
    f_ε =  f_ε[ind]
    ε = ε[ind]
    # get exponent
    @. model(x, p) = p[1]*x+p[2]
    fit = curve_fit(model, vec(log10.(ε)), vec(log10.(f_ε)), [2., 2.])
    D = coef(fit)
    Dϵ = estimate_errors(fit)
    return D[1], Dϵ[1], vec(log10.(ε)),vec(log10.(f_ε))
end



function wel_var(M₂, μ, xₙ, n)

    μ₂ = μ + (xₙ - μ)/n
    M₂ = M₂ + (xₙ - μ)*(xₙ - μ₂)
    return M₂
end

function wel_mean(μ, xₙ, n)
    return μ + (xₙ - μ)/n
end


function static_estimate(xg,yg,bsn)

    xg = bsn.xg; yg = bsn.yg;
    nx=length(xg)
    ny=length(yg)

    N_u = [] # number of uncertain box
    N = [] # number of boxes
    ε = [] # resolution

    # resolution in log scale
    min_ε = 1;
    max_ε = 10;

    r_ε = min_ε:max_ε

    for (k,eps) in enumerate(r_ε)
        Nb=0; Nu=0; μ=0; σ²=0; M₂=0;
        completed = 0;
        # Find uncertain boxes
        while completed == 0
            k1 = rand(1:nx)
            k2 = rand(eps+1:ny-eps)


            c1 = bsn.basin[k1,k2]
            c2 = bsn.basin[k1,k2+eps]
            c3 = bsn.basin[k1,k2-eps]

            if length(unique([c1,Int(c2),Int(c3)]))>1
                Nu = Nu + 1
            end
            Nb += 1

            # Welford's online average estimation and variance of the estimator
            M₂ = wel_var(M₂, μ, Nu/Nb, Nb)
            μ = wel_mean(μ, Nu/Nb, Nb)
            σ² = M₂/Nb

            # Stopping criterion: variance of the estimator of the mean bellow 1e-3
            if Nu > 50 && σ² < 1e-4
                completed = 1
                #@show Nu,Nb,σ²
            end

        end
        push!(N_u,Nu)
        push!(N,Nb)
        push!(ε,eps)
    end

    # uncertain function
    f_ε = N_u./N

    # remove zeros in case there are:
    ind = f_ε .> 0.
    f_ε =  f_ε[ind]
    ε = ε[ind]
    # get exponent
    @. model(x, p) = p[1]*x+p[2]
    fit = curve_fit(model, vec(log10.(ε)), vec(log10.(f_ε)), [2., 2.])
    D = coef(fit)
    #@show estimate_errors(fit)
    return D[1],ε,f_ε

end
