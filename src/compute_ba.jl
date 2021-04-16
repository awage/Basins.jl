

mutable struct basin_info{I,F,V,S}
    basin :: I
    xg :: F
    yg :: F
    iter_f! :: Function
    reinit_f! :: Function
    idxs :: S
    current_color :: Int64
    next_avail_color :: Int64
    consecutive_match :: Int64
    consecutive_other_basins :: Int64
    prevConsecutives :: Int64
    prev_attr :: Int64
    prev_bas :: Int64
    prev_step :: Int64
    step :: Int64
    attractors :: V
end


function get_box(u, bsn_nfo::basin_info)
    xg = bsn_nfo.xg; yg = bsn_nfo.yg; # aliases
    xstep = (xg[2]-xg[1])
    ystep = (yg[2]-yg[1])

    xu=u[1]
    yu=u[2]
    n=0; m=0;
    # check boundary
    if xu >= xg[1] && xu <= xg[end] && yu >= yg[1] && yu <= yg[end]
        n = Int(round((xu-xg[1])/xstep))+1
        m = Int(round((yu-yg[1])/ystep))+1 # +1 for 1 based indexing
    else
        n=-1
        m=-1
    end
    return n,m
end

## Procedure described in  H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012
# The idea is to color the grid with the current color. When an attractor box is hit (even color), the initial condition is colored
# with the color of its basin (odd color). If the trajectory hits another basin 10 times in row the IC is colored with the same
# color as this basin.
function procedure!(bsn_nfo::basin_info, n, m, u)
    max_check = 60
    next_c = bsn_nfo.basin[n,m]
    bsn_nfo.step += 1

    if iseven(next_c) && bsn_nfo.consecutive_match < max_check
        # check wether or not we hit an attractor (even color). Make sure we hit two consecutive times.
        if bsn_nfo.prev_attr == next_c
            bsn_nfo.prevConsecutives +=1
        else
            bsn_nfo.prev_attr = next_c
            bsn_nfo.prevConsecutives =1
            return 0;
        end

        if bsn_nfo.prevConsecutives >= 2
        # Wait if we hit the attractor a second time just to check if it is not a nearby trajectory
            c3 = next_c+1
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            [ bsn_nfo.basin[k[1],k[2]] = c3  for k in ind]
            reset_bsn_nfo!(bsn_nfo)
            return 1
         end
    end

    if next_c == 1 && bsn_nfo.consecutive_match < max_check
        # uncolored box, color it with current odd color
        bsn_nfo.basin[n,m] = bsn_nfo.current_color + 1
        bsn_nfo.consecutive_match = 0
        return 0
    elseif next_c == 1 && bsn_nfo.consecutive_match >= max_check
        # Maybe chaotic attractor, perodic or long recursion.
        bsn_nfo.basin[n,m] = bsn_nfo.current_color
        #println("1 y > max_check")
        return 0
    elseif next_c == bsn_nfo.current_color + 1
        # hit a previously visited box with the current color, possible attractor?
        if bsn_nfo.consecutive_match < max_check
            bsn_nfo.consecutive_match += 1
            return 0
        else
            # println("got attractor")
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            [ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
            bsn_nfo.basin[n,m] = bsn_nfo.current_color
            push!(bsn_nfo.attractors, [bsn_nfo.current_color/2, u]) # store attractor
            # We continue iterating until we hit again the same attractor. In which case we stop.
            return 0;
        end
    elseif isodd(next_c) && 0 < next_c < bsn_nfo.current_color &&  bsn_nfo.consecutive_match < max_check
        # hit a colored basin point of the wrong basin, happens all the time, we check if it happens 10 times in a row or if it happens N times along the trajectory whether to decide if it is another basin.
        bsn_nfo.consecutive_other_basins += 1

        if bsn_nfo.prev_bas == next_c &&  bsn_nfo.prev_step == bsn_nfo.step-1
            bsn_nfo.prevConsecutives +=1
            bsn_nfo.prev_step += 1
        else
            bsn_nfo.prev_bas = next_c
            bsn_nfo.prev_step = bsn_nfo.step
            bsn_nfo.prevConsecutives =1
        end

        if bsn_nfo.consecutive_other_basins > 60 || bsn_nfo.prevConsecutives > 10
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            [bsn_nfo.basin[k[1],k[2]] = next_c for k in ind]
            reset_bsn_nfo!(bsn_nfo)
            return 1
        end
        return 0

    elseif iseven(next_c) && bsn_nfo.consecutive_match >= max_check
        # We have checked the presence of an attractor: tidy up everything and get a new box.
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
        bsn_nfo.current_color = bsn_nfo.next_avail_color
        bsn_nfo.next_avail_color += 2
        #println("even and > max check ", next_c)
        reset_bsn_nfo!(bsn_nfo)
        return 1;
    else
        return 0
    end
end


function check_outside_the_screen!(bsn_nfo::basin_info, new_u, old_u, ni, mi, inlimbo)

    if norm(new_u-old_u) < 1e-5
        #println("Got stuck somewhere, Maybe an attractor outside the screen: ", new_u)
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        reset_bsn_nfo!(bsn_nfo)
        bsn_nfo.basin[ni,mi]=-1 # this CI goes to a attractor outside the screen, set to -1 (even color)
        return -1  # get next box
    elseif inlimbo > 60*20
        #println("trajectory diverges: ", new_u)
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        reset_bsn_nfo!(bsn_nfo)
        bsn_nfo.basin[ni,mi]=-1 # this CI is problematic or diverges, set to -1 (even color)
        return -1  # get next box
    end
    return 0
end

# This routine checks what the is doing the attractor outside the screen, do not modify basin.
function check_outside_the_screen(new_u, old_u, inlimbo)

    if norm(new_u-old_u) < 1e-5
        #println("Got stuck somewhere, Maybe an attractor outside the screen: ", new_u)
        # this CI goes to a attractor outside the screen, set to -1 (even color)
        return -1  # get next box
    elseif inlimbo > 60*20
        #println("trajectory diverges: ", new_u)
        # this CI is problematic or diverges, set to -1 (even color)
        return -1  # get next box
    end
    return 0
end


function reset_bsn_nfo!(bsn_nfo::basin_info)
    #@show bsn_nfo.step
    bsn_nfo.consecutive_match = 0
    bsn_nfo.consecutive_other_basins = 0
    bsn_nfo.prevConsecutives = 0
    bsn_nfo.prev_attr = 1
    bsn_nfo.prev_bas = 1
    bsn_nfo.prev_step = 0
    bsn_nfo.step = 0
end

"""
    draw_basin(xg, yg, integ, iter_f!::Function, reinit_f!::Function)
Compute an estimate of the basin of attraction on a two-dimensional plane. Low level function,
for higher level functions see: `basin_poincare_map`, `basin_discrete_map`, `basin_stroboscopic_map`

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.
* `iter_f!` : function that iterates the map or the system, see step! from DifferentialEquations.jl and
examples for a Poincaré map of a continuous system.
* `reinit_f!` : function that sets the initial condition to test.
"""
function draw_basin(xg, yg, integ, iter_f!::Function, reinit_f!::Function; idxs=1:2)


    complete = 0;

    bsn_nfo = basin_info(ones(Int16, length(xg), length(yg)), xg, yg, iter_f!, reinit_f!, idxs, 2,4,0,0,0,1,1,0,0,[])

    reset_bsn_nfo!(bsn_nfo)

    while complete == 0
         # pick the first empty box
         get_empt_box = findall(bsn_nfo.basin .== 1)
         if length(get_empt_box) == 0
             complete = 1
             break
         end

         ni = get_empt_box[1][1]
         mi = get_empt_box[1][2]
         x0 = xg[ni]
         y0 = yg[mi]

         # Tentatively assign a color: odd is for basins, even for attractors.
         # First color is one
         bsn_nfo.basin[ni,mi] = bsn_nfo.current_color + 1

         # reinitialize integrator
         u0 = [x0, y0]
         reinit_f!(integ,u0)
         next_box = 0
         inlimbo = 0

         while next_box == 0
            old_u = integ.u[idxs]
            iter_f!(integ) # perform a step
            new_u = integ.u[idxs]
            n,m = get_box(new_u, bsn_nfo)
            if n>=0 # apply procedure only for boxes in the defined space
                next_box = procedure!(bsn_nfo, n, m, new_u)
                inlimbo = 0
            else
                # We are outside the box: anything can happen!
                inlimbo +=1
            end

            if inlimbo > 60
                # If we stay too long outside we check if the trajectory is stuck at some unknown attractor outside.
                next_box = check_outside_the_screen!(bsn_nfo, new_u, old_u, ni, mi, inlimbo)
            end

         end
    end


    ind  = findall(iseven.(bsn_nfo.basin) .== true)
    [bsn_nfo.basin[k] = bsn_nfo.basin[k]+1 for k in ind ]
    bsn_nfo.basin = (bsn_nfo.basin .- 1).//2

    return bsn_nfo
end


"""
    basin_poincare_map(xg, yg, integ; kwargs...)
Compute an estimate of the basin of attraction on a two-dimensional plane using a Poincaré map.

[H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012]

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.

## Keyword Arguments
* `plane` A `Tuple{Int, <: Number}`, like `(j, r)` : the plane is defined
  as when the `j` variable of the system equals the value `r`. It can also be
  a vector of length `D+1`. The first `D` elements of the
  vector correspond to ``\\mathbf{a}`` while the last element is ``b``. See ChaosTools.jl
* `Tmax` : maximum time to search for an intersection with the plane before giving up.
* `direction = -1` : Only crossings with `sign(direction)` are considered to be long to the surface of section.
Positive direction means going from less than ``b`` to greater than ``b``.
* `idxs = 1:D` : Optionally you can choose which variables to save. Defaults to
the entire state.
* `rootkw = (xrtol = 1e-6, atol = 1e-6)` : A `NamedTuple` of keyword arguments
passed to `find_zero` from [Roots.jl](https://github.com/JuliaMath/Roots.jl).
"""
function basin_poincare_map(xg, yg, integ; plane=(3,0.), Tmax = 20.,
    direction = -1, idxs = 1:2, rootkw = (xrtol = 1e-6, atol = 1e-6))

    if length(idxs) > 2
        @error "Can only compute basins in two dimensions"
    end

    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)
    planecrossing = PlaneCrossing(plane, direction > 0)
    f = (t) -> planecrossing(integ(t))

    # set the iterator function with the low level function for the Poincaré Map
    iter_f! = (integ) -> DynamicalSystems.ChaosTools.poincaremap!(integ, f, planecrossing, integ.t+Tmax, i, rootkw)

    # Carefully set the initial conditions on the defined plane and
    reinit_f! = (integ,y) -> _initf(integ, y, idxs, plane)

    basin = draw_basin(xg, yg, integ, iter_f!, reinit_f!; idxs=i)


end

function _initf(integ, y, idxs, plane)
    j = plane[1]; v = plane[2]; # take care of the plane
    u = zeros(1,length(integ.u))
    u[idxs] = y
    u[j] = v
    # all other coordinates are zero
    reinit!(integ, u)
end





"""
    basin_stroboscopic_map(xg, yg, integ; T=1., idxs=1:2)
Compute an estimate of the basin of attraction on a two-dimensional plane using a stroboscopic map.

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.

## Keyword Arguments
* `T` : Period of the stroboscopic map
* `idxs = 1:D` : Optionally you can choose which variables to save. Defaults to the entire state.
"""
function basin_stroboscopic_map(xg, yg, integ; T=1., idxs=1:2)

    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)

    iter_f! = (integ) -> step!(integ, T, true)
    reinit_f! =  (integ,y) -> _init(integ, y, i)

    return draw_basin(xg, yg, integ, iter_f!, reinit_f!; idxs=i)
end

function _init(integ, y, idxs)
    u = zeros(length(integ.u))
    u[idxs] = y
    # all other coordinates are zero
    reinit!(integ, u)
end

"""
    basin_discrete_map(xg, yg, integ; idxs=1:2)
Compute an estimate of the basin of attraction on a two-dimensional plane using a discrete map.

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator of the discrete system

## Keyword Arguments
* `idxs = 1:D` : Optionally you can choose which variables to save. Defaults to the entire state.
"""
function basin_discrete_map(xg, yg, integ; idxs=1:2)

    iter_f! = (integ) -> step!(integ)
    reinit_f! =  (integ,y) -> _init(integ, y, idxs)
    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)

    return draw_basin(xg, yg, integ, iter_f!, reinit_f!; idxs=i)
end



# identify an attractor from the basin.
function get_IC_color!(bsn_nfo::basin_info, n,m)
    bsn_nfo.step += 1
    next_c = bsn_nfo.basin[n,m]

        # This routine aims at iditenfying an attractor, when the same bassin has been hit 10 times in a row we assume the original IC belong to this bassin.

        if bsn_nfo.prev_bas == next_c &&  bsn_nfo.prev_step == bsn_nfo.step-1
            bsn_nfo.prevConsecutives +=1
            bsn_nfo.prev_step += 1
        else
            bsn_nfo.prev_bas = next_c
            bsn_nfo.prev_step = bsn_nfo.step
            bsn_nfo.prevConsecutives =1
        end

        if bsn_nfo.prevConsecutives >=10
            reset_bsn_nfo!(bsn_nfo)
            return next_c
        end
        return 0
end


function get_color_point!(bsn_nfo::basin_info, integ, u0)
    # This routine identifies the attractor using the previously defined basin.
    # The idea is that is we hit the same basin 10 times

    # reinitialize integrator
    bsn_nfo.reinit_f!(integ, u0)
    reset_bsn_nfo!(bsn_nfo)

    done = 0;
    inlimbo = 0

    while done == 0
       old_u = integ.u[bsn_nfo.idxs]
       bsn_nfo.iter_f!(integ)
       new_u = integ.u[bsn_nfo.idxs]

       n,m = get_box(new_u, bsn_nfo)

       if n>=0 # apply procedure only for boxes in the defined space
           done = get_IC_color!(bsn_nfo, n, m)
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



# Follow the trajectory until we hit the attractor
function get_color_precise!(bsn_nfo::basin_info, integ, u0; radius=0.)

    # reinitialize integrator
    bsn_nfo.reinit_f!(integ, u0)
    reset_bsn_nfo!(bsn_nfo)

    # set the radius to the grid resolution. This is a conservative distance
    radius = (radius == 0.) ?  maximum([bsn_nfo.xg[2]-bsn_nfo.xg[1] , bsn_nfo.yg[2]-bsn_nfo.yg[1]]) : radius

    # function for attractor metrics:
    get_min_dist = (u) -> minimum([norm(u-p[2]) for p in bsn_nfo.attractors])

    it_cnt = 0;
    done = 0;
    inlimbo = 0

    while done == 0
       old_u = integ.u[bsn_nfo.idxs]
       integ.t = 0
       bsn_nfo.iter_f!(integ)
       new_u = integ.u[bsn_nfo.idxs]
       n,m = get_box(new_u, bsn_nfo)

       if n>=0 # apply procedure only for boxes in the defined space
           #@show get_min_dist(new_u)
           if get_min_dist(new_u) < radius
               # find the attractor:
               v = [norm(new_u - att[2]) for att in bsn_nfo.attractors]
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

       it_cnt += 1

       if it_cnt > 5000
           @warn "Max iteration in get_color_precise! something went wrong, check the location of attractors."
           break;
       end
    end

    return done
end



"""
    compute_basin_precise(bsn_nfo::basin_info, integ)
Compute the basin of attraction using only the attractors. The attractors must have been found previously in `bsn_nfo.attractors`

## Arguments
* `integ` : the matrix containing the information of the basin.
* `bsn_nfo` : structure that holds the information of the basin as well as the map function. This structure is set when the basin is first computed with `basin_stroboscopic_map` or `basin_poincare_map`.

## Keyword arguments:
* `radius` : this the distance to the attractor in order to identify correctly the initial condition. Default value is the grid size.
"""
function compute_basin_precise(bsn_nfo::basin_info, integ; radius=0.)
    # Experiment: compute the basin when the attractors are known
    new_bas = zeros(Int16, size(bsn_nfo.basin))

    for (k,x) in enumerate(bsn_nfo.xg), (n,y) in enumerate(bsn_nfo.yg)
        new_bas[k,n]=get_color_precise!(bsn_nfo, integ,[x,y]; radius=radius)
    end

    return new_bas
end









## Procedure described in  H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012
# The idea is to color the grid with the current color. When an attractor box is hit (even color), the initial condition is colored
# with the color of its basin (odd color). If the trajectory hits another basin 10 times in row the IC is colored with the same
# color as this basin.
function procedure2!(bsn_nfo::basin_info, n, m, u)
    max_check = 60
    next_c = bsn_nfo.basin[n,m]
    bsn_nfo.step += 1

    if iseven(next_c) && bsn_nfo.consecutive_match < max_check
        # check wether or not we hit an attractor (even color). Make sure we hit two consecutive times.
        if bsn_nfo.prev_attr == next_c
            bsn_nfo.prevConsecutives +=1
        else
            bsn_nfo.prev_attr = next_c
            bsn_nfo.prevConsecutives =1
            return 0;
        end

        if bsn_nfo.prevConsecutives >= 10
        # Wait if we hit the attractor a 10 times in a row just to check if it is not a nearby trajectory
            c3 = next_c+1
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
            reset_bsn_nfo!(bsn_nfo)
            return c3
         end
    end

    if next_c == 1 && bsn_nfo.consecutive_match < max_check
        # uncolored box, color it with current odd color
        bsn_nfo.basin[n,m] = bsn_nfo.current_color + 1
        bsn_nfo.consecutive_match = 0
        return 0
    elseif next_c == 1 && bsn_nfo.consecutive_match >= max_check
        # Maybe chaotic attractor, perodic or long recursion.
        bsn_nfo.basin[n,m] = bsn_nfo.current_color
        #println("1 y > max_check")
        return 0
    elseif next_c == bsn_nfo.current_color + 1
        # hit a previously visited box with the current color, possible attractor?
        if bsn_nfo.consecutive_match < max_check
            bsn_nfo.consecutive_match += 1
            return 0
        else
            println("got attractor")
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            [ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
            bsn_nfo.basin[n,m] = bsn_nfo.current_color
            push!(bsn_nfo.attractors, [bsn_nfo.current_color/2, u]) # store attractor
            # We continue iterating until we hit again the same attractor. In which case we stop.
            return 0;
        end
    elseif iseven(next_c) && bsn_nfo.consecutive_match >= max_check
        # We have checked the presence of an attractor: tidy up everything and get a new box.
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
        bsn_nfo.current_color = bsn_nfo.next_avail_color
        bsn_nfo.next_avail_color += 2
        #println("even and > max check ", next_c)
        reset_bsn_nfo!(bsn_nfo)
        return 1;
    else
        return 0
    end
end


"""
    draw_basin2(xg, yg, integ, iter_f!::Function, reinit_f!::Function)
Compute an estimate of the basin of attraction on a two-dimensional plane. Low level function,
for higher level functions see: `basin_poincare_map`, `basin_discrete_map`, `basin_stroboscopic_map`

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.
* `iter_f!` : function that iterates the map or the system, see step! from DifferentialEquations.jl and
examples for a Poincaré map of a continuous system.
* `reinit_f!` : function that sets the initial condition to test.
"""
function draw_basin2(xg, yg, integ; dt=0.1, idxs=1:2)


    complete = 0;

    reinit_f! = function reinit_f!(integ,u0)
        v0 = zeros(1,dim)
        v0[idxs]=u0
        reinit!(integ,v0)
    end

    iter_f! = (integ) -> step!(integ, T, true)

    bsn_nfo = basin_info(ones(Int16, length(xg), length(yg)), xg, yg, iter_f!, reinit_f!, idxs, 2,4,0,0,0,1,1,0,0,[])
    dim = length(integ.u)
    reset_bsn_nfo!(bsn_nfo)

    while complete == 0
         # pick the first empty box
         get_empt_box = findall(bsn_nfo.basin .== 1)
         if length(get_empt_box) == 0
             complete = 1
             break
         end

         ni = get_empt_box[1][1]
         mi = get_empt_box[1][2]
         x0 = xg[ni]
         y0 = yg[mi]

         # Tentatively assign a color: odd is for basins, even for attractors.
         # First color is one
         bsn_nfo.basin[ni,mi] = bsn_nfo.current_color + 1

         # reinitialize integrator
         u0 = zeros(1,dim)
         u0[idxs]=[x0, y0]
         reinit!(integ,u0)
         next_box = 0
         inlimbo = 0

         while next_box == 0
            old_u = integ.u[idxs]
            step!(integ, dt, true) # perform a step
            new_u = integ.u[idxs]
            n,m = get_box(new_u, bsn_nfo)
            if n>=0 # apply procedure only for boxes in the defined space
                next_box = procedure2!(bsn_nfo, n, m, new_u)
                if next_box > 1
                    bsn_nfo.basin[ni,mi] = next_box
                end
                inlimbo = 0
            else
                # We are outside the box: anything can happen!
                inlimbo +=1
            end

            if inlimbo > 60
                # If we stay too long outside we check if the trajectory is stuck at some unknown attractor outside.
                next_box = check_outside_the_screen!(bsn_nfo, new_u, old_u, ni, mi, inlimbo)
            end

         end
    end

    ind  = findall(iseven.(bsn_nfo.basin) .== true)
    [bsn_nfo.basin[k] = bsn_nfo.basin[k]+1 for k in ind ]
    bsn_nfo.basin = (bsn_nfo.basin .- 1).//2

    return bsn_nfo
end
