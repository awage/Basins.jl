



function bisection_refine!(u_A, u_B, bsn_nfo, integ, basin_A, basin_B, tol)

    # shortcut functions
    function get_col(u0)
        a = get_color_point!(bsn_nfo, integ, u0)
        return iseven(a) ? Int(a/2) : Int((a-1)/2)
    end

    function get_col_next(u0)
        bsn_nfo.reinit_f!(integ,u0)
        bsn_nfo.iter_f!(integ)
        a=get_col(bsn_nfo.get_u(integ))
        return  a
    end

    # directing vector
    Δu = u_B - u_A

    # Midpoint method
    u_A_next=u_A
    u_B_next=u_B
    dist= norm(u_A_next-u_B_next)

    if dist < tol

        if get_col_next(u_A_next) ∉ basin_A
            # This basin has switch we try to perturb u_A_next until we get back to the initial basin.
            # We guarantee that the next iteration will stay also in the good basin
            eps=1e-8
            u_A_rnd = u_A_next + eps*randn(2)
            while get_col_next(u_A_rnd) ∉ basin_A
                eps += eps
                @show u_A_rnd = u_A_next + eps*randn(2)
            end
            u_A_next=u_A_rnd
        end

        if get_col_next(u_B_next) ∉ basin_B
            # This basin has switch we try to perturb u_B_next until we get back to the initial basin.
            eps=1e-8
            u_B_rnd = u_B_next + eps*randn(2)
            while get_col_next(u_B_rnd) ∉ basin_B
                eps += eps
                @show u_B_rnd = u_B_next + eps*randn(2)
            end
            u_B_next=u_B_rnd
        end

        return u_A_next,u_B_next

    else
        # refinement with midpoint method.
        while dist > tol

            mV= Δu/2

            if get_col(u_A_next+mV) ∈ basin_A
                u_A_next = u_A_next+mV
            elseif get_col(u_B_next-mV) ∈ basin_B
                u_B_next = u_B_next-mV
            else
                @show get_col(u_B_next-mV), get_col(u_A_next+mV)
                break
            end

            if get_col_next(u_A_next) ∉ basin_A
                # This basin has switch we try to perturb u_A_next until we get back to the initial basin.
                eps=1e-8
                u_A_rnd = u_A_next + eps*randn(2)
                while get_col_next(u_A_rnd) ∉ basin_A
                    eps += eps
                    @show u_A_rnd = u_A_next + eps*randn(2)
                end
                u_A_next=u_A_rnd
            end


            if get_col_next(u_B_next) ∉ basin_B
                # This basin has switch we try to perturb u_B_next until we get back to the initial basin.
                eps=1e-8
                u_B_rnd = u_B_next + eps*randn(2)
                while get_col_next(u_B_rnd) ∉ basin_B
                    eps += eps
                    @show u_B_rnd = u_B_next + eps*randn(2)
                end
                u_B_next=u_B_rnd
            end

            Δu = u_B_next - u_A_next
            dist= norm(u_A_next-u_B_next)
        end
    end

    return u_A_next,u_B_next
end



"""
    compute_saddle(bsn_nfo::BasinInfo, integ; max_iter=10)
The algorithm the saddle that lies in  a boundary of the basin of attraction of a dynamical system. When the boundary is fractal, this
set is known as the chaotic saddle and is responsible for the transient dynamics of the system. The saddle is computed with the Saddle-Straddle
algorithm. It is necessary to define two `generalized basin`, that is we must separate the basin into two sets (see also keyword arguments).

[H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012]

## Arguments
* `bsn_nfo` : structure that holds the information of the basin as well as the map function. This structure is set when the basin is first computed with `basin_stroboscopic_map` or `basin_poincare_map`.
* `integ` : the matrix containing the information of the basin.
* `bas_A` : vector with the indices of the attractors that will represent the generalized basin A
* `bas_B` : vector with the indices of the attractors that will represent the generalized basin B. Notice that `bas_A ∪ bas_B = [1:N]` and `bas_A ∩ bas_B = ∅`

## Keyword arguments
* `N` : number of points of the saddle to compute
* `init_tol`: length of the straddle segment

## Example
```
# compute 1000 points between generalized A = [1] and basin B = [2,3]
sa,sb = compute_saddle(bsn, integ_df, [1], [2,3],1000)
# sa is the "left" set and sb the right "set"
```

"""
function compute_saddle(integ, bsn_nfo::BasinInfo, bas_A, bas_B; N=100, init_tol = 1e-6)

    # shortcut functions
    function get_col(u0)
        a = get_color_point!(bsn_nfo, integ, u0)
        return iseven(a) ? Int(a/2) : Int((a-1)/2)
    end

    Na = bsn_nfo.Na
    # basic check
    if (Set(bas_A) ∪ Set(bas_B)) != Set(collect(1:Na))
        @error "Generalized basins are not well defined"
        return nothing,nothing
    end

    if length(Set(bas_A) ∩ Set(bas_B)) > 0
        @error "Generalized basins are not well defined"
        return nothing, nothing
    end


    tol = init_tol
    xg = bsn_nfo.xg; yg = bsn_nfo.yg; # aliases

    # Set initial condition near attractors of the selected basins.
    u_A = u_B = [0., 0.]
    v = bsn_nfo.attractors[bas_A[1]*2]
    u_A = v[1]
    v = bsn_nfo.attractors[bas_B[1]*2]
    u_B = v[1]


    saddle_series_A = Array{typeof(u_A),1}(undef, N)
    saddle_series_B = Array{typeof(u_A),1}(undef, N)
    u_A_it = u_B_it = [0.,0.]
    Ttr= 20 # skip 20 points
    switch_bas = 0 # keep track of switching between basin if it happens too much the basin is problematic
    k=1;
    while k<= N

        # Initial bisection
        u_A_r,u_B_r = bisection_refine!(u_A, u_B, bsn_nfo, integ, bas_A, bas_B, tol)

        dist = norm(u_A_r-u_B_r)
        bsn_nfo.reinit_f!(integ,u_A_r)
        bsn_nfo.iter_f!(integ)
        u_A_it = deepcopy(bsn_nfo.get_u(integ)) # for some reason I have to deepcopy this vector
        ca=get_col(u_A_it)

        bsn_nfo.reinit_f!(integ,u_B_r)
        bsn_nfo.iter_f!(integ)
        u_B_it = deepcopy(bsn_nfo.get_u(integ))
        cb=get_col(u_B_it)


        if ca ∉ bas_A || cb ∉ bas_B
            # one of the basins switched attractor, we retry refinement but with bigger tolerance.
            println("Basin have switched this is unusual, basin A:", ca, " basin B:  ", cb)
            switch_bas += 1
            tol = tol*2
            Ttr= 20; # skip 20 points
        else
            if Ttr <= 0 &&  dist < tol
                saddle_series_A[k]=u_A
                saddle_series_B[k]=u_B
                k +=1
                switch_bas = 0
            else
                Ttr -= 1;
            end
            u_A=u_A_it;
            u_B=u_B_it;
            #@show u_A,u_B
            tol = init_tol
        end

        if switch_bas > 20
            @warn "Problem with the basin: it may be a riddled boundary or the initial resolution of the basin is too low."
            break
        end
    end

    return saddle_series_A, saddle_series_B
end
