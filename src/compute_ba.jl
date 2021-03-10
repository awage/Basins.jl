


mutable struct basin_info
    basin :: Array{Int16,2}
    xg :: Array{Float64,1}
    yg :: Array{Float64,1}
    current_color :: Int32
    next_avail_color :: Int32
    consecutive_match :: Int32
    consecutive_other_basins :: Int32
    prevConsecutives :: Int32
    prev_attr :: Int32
    prev_bas :: Int32
    prev_step :: Int32
    step :: Int32

end

# Helper function to optimize allocations
function find_and_replace!(basin, old_c, new_c)
    for (k,c) in enumerate(basin)
        if c == old_c
            basin[k] = new_c
        end
    end
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


function procedure!(bsn_nfo::basin_info, n,m)
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
            #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            #[ bsn_nfo.basin[k[1],k[2]] = c3  for k in ind]
            find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, c3)

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
            #println("got attractor")
            #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            #[ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
            find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, 1)

            bsn_nfo.basin[n,m] = bsn_nfo.current_color
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
            #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
            #[bsn_nfo.basin[k[1],k[2]] = next_c for k in ind]
            find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, next_c)

            reset_bsn_nfo!(bsn_nfo)
            return 1
        end
        return 0

    elseif iseven(next_c) && bsn_nfo.consecutive_match >= max_check
        # We have checked the presence of an attractor: tidy up everything and get a new box.
        #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        #[ bsn_nfo.basin[k[1],k[2]] = 1 for k in ind]
        find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, 1)
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
        #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        #[ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, 1)
        reset_bsn_nfo!(bsn_nfo)
        bsn_nfo.basin[ni,mi]=-1 # this CI goes to a attractor outside the screen, set to -1 (even color)
        return -1  # get next box
    elseif inlimbo > 60*20
        #println("trajectory diverges: ", new_u)
        #ind = findall(bsn_nfo.basin .== bsn_nfo.current_color+1)
        #[ bsn_nfo.basin[k[1],k[2]] = 1  for k in ind]
        find_and_replace!(bsn_nfo.basin, bsn_nfo.current_color+1, 1)
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



function draw_basin(xg, yg, integ_df; T=0.01)

    complete = 0;

    bsn_nfo = basin_info(ones(Int8, length(xg), length(yg)), xg, yg, 2,4,0,0,0,1,1,0,0)

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

         # Tentatively assign a color
         bsn_nfo.basin[ni,mi] = bsn_nfo.current_color + 1

         # reinitialize integrator
         u0 = [x0, y0]
         reinit!(integ_df, u0, t0=0, erase_sol=true,reinit_callbacks=true)

         next_box = 0
         inlimbo = 0

         while next_box == 0
            old_u = integ_df.u
            step!(integ_df, T, true) # perform a step
            new_u = integ_df.u
            n,m = get_box(new_u, bsn_nfo)
            if n>=0 # apply procedure only for boxes in the defined space
                next_box = procedure!(bsn_nfo, n, m)
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
    return bsn_nfo.basin
end


# identify an attractor
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

        if bsn_nfo.prevConsecutives >=5
            reset_bsn_nfo!(bsn_nfo)
            return next_c
        end
        return 0
end


function get_color_point!(bsn_nfo::basin_info, integ, x0, y0, T)
    # This routine identifies the attractor using the previously defined basin.
    # The idea is that is we hit the same basin 10 times

    # reinitialize integrator
    u0 =  [x0, y0]
    reinit!(integ, u0, t0=0, erase_sol=true)
    reset_bsn_nfo!(bsn_nfo)

    done = 0;
    inlimbo = 0

    while done == 0
       old_u = integ.u
       step!(integ, T, true)
       new_u = integ.u

       n,m = get_box(new_u, bsn_nfo)

       if n>=0 # apply procedure only for boxes in the defined space
           done = get_IC_color!(bsn_nfo, n, m)
           inlimbo = 0
       else
           inlimbo +=1
       end

       if inlimbo > 60
           done = check_outside_the_screen(new_u, old_u, inlimbo)
       end

    end
    #@show done
    return done
end
