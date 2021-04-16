

mutable struct basin_info{I,F,V}
    basin :: I
    xg :: F
    yg :: F
    iter_f! :: Function
    reinit_f! :: Function
    get_u :: Function
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


"""
    basin_poincare_map(xg, yg, integ; kwargs...)
Compute an estimate of the basin of attraction on a two-dimensional plane using a Poincaré map.

[H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 1997]

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
function basin_poincare_map(xg, yg, pmap, idxs = 1:2; Ncheck = 2)
    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)
    # Carefully set the initial conditions on the defined plane and
    reinit_f! = (integ,y) -> _initf(integ, y, i)
    get_u = (pmap) -> pmap.integ.u[i]
    basin = draw_basin!(xg, yg, pmap, step!, reinit_f!, get_u, Ncheck)

end


function _initf(pmap, y, idxs)
    u = zeros(1,length(pmap.integ.u))
    u[idxs] = y
    # all other coordinates are zero
    reinit!(pmap, u)
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
function basin_stroboscopic_map(xg, yg, integ; T=1., idxs=1:2, Ncheck = 2)

    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)
    iter_f! = (integ) -> step!(integ, T, true)
    reinit_f! =  (integ,y) -> _init(integ, y, i)
    get_u = (integ) -> integ.u[i]
    return draw_basin!(xg, yg, integ, iter_f!, reinit_f!,get_u, Ncheck)

end

function basin_general_ds(xg, yg, integ; dt=1., idxs=1:2, Ncheck = 10)

    return basin_stroboscopic_map(xg, yg, integ; T=dt, idxs=idxs, Ncheck = Ncheck)

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

    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)

    iter_f! = (integ) -> step!(integ)
    reinit_f! =  (integ,y) -> _init(integ, y, i)
    get_u = (integ) -> integ.u[i]

    return draw_basin!(xg, yg, integ, iter_f!, reinit_f!, get_u, 2)
end






## Procedure described in  H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 2012
# The idea is to color the grid with the current color. When an attractor box is hit (even color), the initial condition is colored
# with the color of its basin (odd color). If the trajectory hits another basin 10 times in row the IC is colored with the same
# color as this basin.
function procedure!(bsn_nfo::basin_info, n, m, u, Ncheck)
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

        if bsn_nfo.prevConsecutives >= Ncheck
            # Wait if we hit the attractor a Ncheck times in a row just to check if it is not a nearby trajectory
            #println("found IC")
            c3 = next_c+1
            ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            if Ncheck == 2
                [ bsn_nfo.basin[k[1],k[2]] = c3  for k in ind]
            else
                [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind] # erase visited boxes
            end
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
        bsn_nfo.consecutive_match = max_check
        println("1 y > max_check")
        return 0
    elseif next_c == bsn_nfo.current_color + 1
        # hit a previously visited box with the current color, possible attractor?
        if bsn_nfo.consecutive_match < max_check
            bsn_nfo.consecutive_match += 1
            return 0
        else
            println("got attractor")
            #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            #[ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
            bsn_nfo.basin[n,m] = bsn_nfo.current_color
            push!(bsn_nfo.attractors, [bsn_nfo.current_color/2, u]) # store attractor
            # We continue iterating until we hit again the same attractor. In which case we stop.
            return 0;
        end
    elseif isodd(next_c) && 0 < next_c < bsn_nfo.current_color &&  bsn_nfo.consecutive_match < max_check && Ncheck == 2
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
            return next_c
        end
        return 0
    elseif iseven(next_c) &&   (max_check <= bsn_nfo.consecutive_match < 2*max_check)
        # We make sure we hit the basin 60 consecutive times
        bsn_nfo.consecutive_match+=1
        return 0
    elseif iseven(next_c) && bsn_nfo.consecutive_match >= max_check*2
        # We have checked the presence of an attractor: tidy up everything and get a new box.
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
        bsn_nfo.basin[n,m] = bsn_nfo.current_color
        bsn_nfo.current_color = bsn_nfo.next_avail_color
        bsn_nfo.next_avail_color += 2
        println("even and > max check ", next_c)
        reset_bsn_nfo!(bsn_nfo)
        return next_c+1;
    else
        return 0
    end
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
function draw_basin!(xg, yg, integ, iter_f!::Function, reinit_f!::Function, get_u::Function, Ncheck)

    complete = 0;

    bsn_nfo = basin_info(ones(Int16, length(xg), length(yg)), xg, yg, iter_f!, reinit_f!, get_u, 2,4,0,0,0,1,1,0,0,[])

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

         u0=[x0, y0]

         bsn_nfo.basin[ni,mi] = get_color_point!(bsn_nfo, integ, u0; Ncheck=Ncheck)
    end

    return bsn_nfo
end



function get_color_point!(bsn_nfo::basin_info, integ, u0; Ncheck=2)
    # This routine identifies the attractor using the previously defined basin.
    # reinitialize integrator
    bsn_nfo.reinit_f!(integ, u0)
    reset_bsn_nfo!(bsn_nfo)

    done = 0;
    inlimbo = 0

    while done == 0
       old_u = bsn_nfo.get_u(integ)
       bsn_nfo.iter_f!(integ)
       new_u = bsn_nfo.get_u(integ)

       n,m = get_box(new_u, bsn_nfo)

       if n>=0 # apply procedure only for boxes in the defined space
           done = procedure!(bsn_nfo, n, m, new_u, Ncheck)
           inlimbo = 0
       else
           # We are outside the defined grid
           inlimbo +=1
       end

       if inlimbo > 60
           done = check_outside_the_screen!(bsn_nfo, new_u, old_u, inlimbo)
       end
    end

    return done
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


function check_outside_the_screen!(bsn_nfo::basin_info, new_u, old_u, inlimbo)

    if norm(new_u-old_u) < 1e-5
        #println("Got stuck somewhere, Maybe an attractor outside the screen: ", new_u)
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        reset_bsn_nfo!(bsn_nfo)
        # this CI goes to a attractor outside the screen, set to -1 (even color)
        return -1  # get next box
    elseif inlimbo > 60*20
        #println("trajectory diverges: ", new_u)
        ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        [ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        reset_bsn_nfo!(bsn_nfo)
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
