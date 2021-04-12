

function bisection_refine!(u_A, u_B, bsn_nfo, integ, basin_A, basin_B, tol)

    # shortcut function
    get_col = (u0) -> get_color_precise!(bsn_nfo, integ, u0)

    function get_col_next(u0)
        bsn_nfo.reinit_f!(integ,u0)
        bsn_nfo.iter_f!(integ)
        return get_col(integ.u)
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
    compute_saddle(integ, bsn_nfo::basin_info; max_iter=10)
The algorithm the saddle that lies in  a boundary of the basin of attraction of a dynamical system. When the boundary is fractal, this
set is known as the chaotic saddle and is responsible for the transient dynamics of the system. The saddle is computed with the Saddle-Straddle
algorithm. It is necessary to define two `generalized basin`, that is we must separate the basin into two sets (see also keyword arguments).

[H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012]

## Arguments
* `integ` : the matrix containing the information of the basin.
* `bsn_nfo` : structure that holds the information of the basin as well as the map function. This structure is set when the basin is first computed with `basin_stroboscopic_map` or `basin_poincare_map`.
* `bas_A` : vector with the indices of the attractors that will represent the generalized basin A
* `bas_B` : vector with the indices of the attractors that will represent the generalized basin B. Notice that `bas_A ∪ bas_B = [1:N]` and `bas_A ∩ bas_B = ∅`

## Keyword arguments
* `N` : number of points of the saddle to compute

## Example
```
# compute 1000 points between generalized A = [1] and basin B = [2,3]
sa,sb = compute_saddle(bsn, integ_df, [1], [2,3],1000)
# sa is the "left" set and sb the right "set"
```

"""
function compute_saddle(bsn_nfo::basin_info, integ, bas_A, bas_B; N=100)


    x0=0.
    tol = init_tol = 1e-6

    xg = bsn_nfo.xg; yg = bsn_nfo.yg; # aliases

    # Set initial condition near attractors of the selected basins.
    u_A = u_B = [0., 0.]
    for p in bsn_nfo.attractors
        @show  p
        if p[1] == bas_A[1]
            u_A = [p[2], p[3]]
        end
        if p[1] == bas_B[1]
            u_B = [p[2], p[3]]
        end
    end


    @show get_color_precise!(bsn_nfo,integ,u_A) ∈ bas_A

    @show get_color_precise!(bsn_nfo,integ,u_B) ∈ bas_B

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
        u_A_it = deepcopy(get_state(integ)) # for some reason I have to deepcopy this vector
        ca=get_color_precise!(bsn_nfo,integ,u_A_it)

        bsn_nfo.reinit_f!(integ,u_B_r)
        bsn_nfo.iter_f!(integ)
        u_B_it = deepcopy(get_state(integ))
        cb=get_color_precise!(bsn_nfo,integ,u_B_it)


        if ca ∉ bas_A || cb ∉ bas_B
            # one of the basins switched attractor, we retry refinement but with bigger tolerance.
            println("Basin have switched this is unusual, basin A:", ca, " basin B:  ", cb)
            switch_bas += 1
            tol = tol*2
            Ttr= 20; # skip 20 points
        else
            if Ttr <= 0 &&  dist < tol
                #push!(saddle_series_A, u_A);
                #push!(saddle_series_B, u_B);
                saddle_series_A[k]=u_A
                saddle_series_B[k]=u_B
                k +=1
                switch_bas = 0
            else
                Ttr -= 1;
            end
            u_A=u_A_it;
            u_B=u_B_it;
            tol = init_tol
        end

        if switch_bas > 20
            @warn "Problem with the basin: it may be a riddled boundary or the initial resolution of the basin is too low."
            break
        end
    end

    return saddle_series_A, saddle_series_B
end



# Follow the trajectory until we hit the attractor
function get_color_precise!(bsn_nfo::basin_info, integ, u0)


    # reinitialize integrator
    bsn_nfo.reinit_f!(integ, u0)
    reset_bsn_nfo!(bsn_nfo)
    radius = maximum([bsn_nfo.xg[2]-bsn_nfo.xg[1] , bsn_nfo.yg[2]-bsn_nfo.yg[1]])/2

    # function for attractor metrics:
    get_min_dist = (u) -> minimum([norm(u-[p[2],p[3]]) for p in bsn_nfo.attractors])

    done = 0;
    inlimbo = 0

    while done == 0
       old_u = integ.u
       integ.t = 0
       bsn_nfo.iter_f!(integ)
       new_u = integ.u

       n,m = get_box(new_u, bsn_nfo)

       if n>=0 # apply procedure only for boxes in the defined space
           if get_min_dist(new_u) < radius
               # find the attractor:
               v = [norm(new_u-[p[2],p[3]]) for p in bsn_nfo.attractors]
               for (k,m) in enumerate(v)
                   if get_min_dist(new_u) == m
                       return bsn_nfo.attractors[k][1]
                   end
               end

           end
           inlimbo = 0
       else
           # We are outside the defined grid
           inlimbo +=1
       end

       if inlimbo > 60
           done = check_outside_the_screen(new_u, old_u, inlimbo)
       end

    end

    return done
end

# Experiment: compute the basin when the attractors are known
function compute_basin_precise(bsn_nfo::basin_info, integ)


new_bas = zeros(Int16, size(bsn_nfo.basin))

for (k,x) in enumerate(bsn_nfo.xg), (n,y) in enumerate(bsn_nfo.yg)
    new_bas[k,n]=get_color_precise!(bsn_nfo,integ,[x,y])
end

return new_bas
end
