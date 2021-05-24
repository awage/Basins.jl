PoincareMap = ChaosTools.PoincareMap

mutable struct BasinInfo{I,F,D,Q}
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
    attractors :: D #Dict{Int16,Vector{Vector{Float64}}}
    Na :: Int64
    visited :: Q
    available :: Q
end

function Base.show(io::IO, bsn::BasinInfo)
    println(io, "Basin of attraction structure")
    println(io,  rpad(" size : ", 14),    size(bsn.basin))
    println(io,  rpad(" Number of attractors found: ", 14),   bsn.Na  )
end


"""
    basins_map2D(xg, yg, integ; kwargs...)
Compute an estimate of the basin of attraction on a two-dimensional plane using a map of the plane onto itself.
The dynamical system should be a discrete two dimensional system such as:
    * Discrete 2D map.
    * 2D poincaré map.
    * A 2D stroboscopic map.

[H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 1997]

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : A  integrator handle of the dynamical system. For a Poincaré map, the handle is a `pmap`
as defined in [`poincaremap`](@ref)

## Keyword Arguments
* `T` : Period of the stroboscopic map.
* `Ncheck` : A parameter that sets the number of consecutives hits of an attractor before deciding the basin
of the initial condition.


## Example:

```jl
using DynamicalSystems, ChaosTools
ds = Systems.rikitake(μ = 0.47, α = 1.0)
xg=range(-6.,6.,length=150); yg=range(-6.,6.,length=150)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6; idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9)
bsn = basins_map2D(xg, yg, pmap)
```

"""
function basins_map2D(xg, yg, pmap::PoincareMap; Ncheck = 3)
    reinit_f! = (pmap,y) -> _init_map(pmap, y, pmap.i)
    get_u = (pmap) -> pmap.integ.u[pmap.i]
    basin = draw_basin!(xg, yg, pmap, step!, reinit_f!, get_u, Ncheck)
end


function _init_map(pmap::PoincareMap, y, idxs)
    u = zeros(1,length(pmap.integ.u))
    u[idxs] = y
    # all other coordinates are zero
    reinit!(pmap, u)
end


function basins_map2D(xg, yg, integ; T=0., Ncheck = 2)
    if T>0
        iter_f! = (integ) -> step!(integ, T, true)
    else
        iter_f! = (integ) -> step!(integ)
    end
    reinit_f! =  (integ,y) -> reinit!(integ, y)
    get_u = (integ) -> integ.u

    return draw_basin!(xg, yg, integ, iter_f!, reinit_f!, get_u, Ncheck)
end






"""
    basins_general(xg, yg, integ; T=1., idxs=1:2)
Compute an estimate of the basin of attraction on a two-dimensional plane using a stroboscopic map.

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.

## Keyword Arguments
* `dt` : Time step of the integrator. It recommended to use values above 1.
* `idxs = 1:D` : Optionally you can choose which variables to save. Defaults to the entire state.
* `Ncheck` : A parameter that sets the number of consecutives hits of an attractor before deciding the basin
of the initial condition.

## Example:
```jl
using DynamicalSystems, ChaosTools
ds = Systems.magnetic_pendulum(γ=1, d=0.2, α=0.2, ω=0.8, N=3)
integ = integrator(ds, u0=[0,0,0,0], reltol=1e-9)
xg=range(-4,4,length=150)
yg=range(-4,4,length=150)
@time bsn = basins_general(xg, yg, integ; dt=1., idxs=1:2)
"""
function basins_general(xg, yg, integ; dt=1., idxs=1:2, Ncheck = 10)
    i = typeof(idxs) <: Int ? i : SVector{length(idxs), Int}(idxs...)
    iter_f! = (integ) -> step!(integ, dt, true)
    reinit_f! =  (integ,y) -> _init_ds(integ, y, i)
    get_u = (integ) -> integ.u[i]
    return draw_basin!(xg, yg, integ, iter_f!, reinit_f!,get_u, Ncheck)
end


function _init_ds(integ, y, idxs)
    u = zeros(length(integ.u))
    u[idxs] = y
    # all other coordinates are zero
    reinit!(integ, u)
end


## Procedure described in  H. E. Nusse and J. A. Yorke, Dynamics: numerical explorations, Springer, New York, 1997
# The idea is to color the grid with the current color. When an attractor box is hit (even color), the initial condition is colored
# with the color of its basin (odd color). If the trajectory hits another basin 10 times in row the IC is colored with the same
# color as this basin.
function procedure!(bsn_nfo::BasinInfo, n, u, Ncheck::Int)
    max_check = 60
    #next_c = bsn_nfo.basin[n]
    next_c=get_data(bsn_nfo.basin, n)
    bsn_nfo.step += 1


    if iseven(next_c) && bsn_nfo.consecutive_match < max_check
        # check wether or not we hit an attractor (even color). Make sure we hit Ncheck consecutive times.
        if bsn_nfo.prev_attr == next_c
            bsn_nfo.prevConsecutives +=1
        else
            # reset prevConsecutives
            bsn_nfo.prev_attr = next_c
            bsn_nfo.prevConsecutives =1
            return 0;
        end

        if bsn_nfo.prevConsecutives >= Ncheck
            # Wait if we hit the attractor a Ncheck times in a row just to check if it is not a nearby trajectory
            #println("found IC")
            c3 = next_c+1

            if Ncheck == 2
                # For maps we can color the previous steps as well. Every point of the trajectory lead
                # to the attractor
                find_and_replace!(bsn_nfo, bsn_nfo.current_color+1, c3)
                #ind = (bsn_nfo.basin .== bsn_nfo.current_color+1)
                #bsn_nfo.basin[ind] .= c3
            else
                # For higher dimensions we erase the past iterations.
                find_and_replace!(bsn_nfo, bsn_nfo.current_color+1, 1)
                #bsn_nfo.basin[ind] .= 1
            end
            reset_bsn_nfo!(bsn_nfo)
            return c3
         end
    end

    if next_c == 1 && bsn_nfo.consecutive_match < max_check
        # uncolored box, color it with current odd color
        #bsn_nfo.basin[n] = bsn_nfo.current_color + 1
        set_data!(bsn_nfo.basin, n, bsn_nfo.current_color+1)
        push!(bsn_nfo.visited,u)
        bsn_nfo.consecutive_match = 0
        return 0
    elseif next_c == 1 && bsn_nfo.consecutive_match >= max_check
        # Maybe chaotic attractor, perodic or long recursion.
        # Color this box as part of an attractor
        #bsn_nfo.basin[n] = bsn_nfo.current_color
        set_data!(bsn_nfo.basin, n, bsn_nfo.current_color)
        store_attractor!(bsn_nfo, u)
        bsn_nfo.consecutive_match = max_check
        #println("1 y > max_check")
        return 0
    elseif next_c == bsn_nfo.current_color + 1
        # hit a previously visited box with the current color, possible attractor?
        if bsn_nfo.consecutive_match < max_check
            bsn_nfo.consecutive_match += 1
            return 0
        else
            #println("got attractor")
            #bsn_nfo.basin[n] = bsn_nfo.current_color
            set_data!(bsn_nfo.basin, n, bsn_nfo.current_color)
            store_attractor!(bsn_nfo, u)
            # We continue iterating until we hit again the same attractor. In which case we stop.
            return 0;
        end
    elseif isodd(next_c) && 0 < next_c < bsn_nfo.current_color &&  bsn_nfo.consecutive_match < max_check && Ncheck == 2
        # hit a colored basin point of the wrong basin, happens all the time, we check if it happens
        #10 times in a row or if it happens N times along the trajectory whether to decide if it is another basin.
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
            #ind = (bsn_nfo.basin .== bsn_nfo.current_color+1)
            #bsn_nfo.basin[ind] .= next_c
            find_and_replace!(bsn_nfo, bsn_nfo.current_color+1, next_c)
            reset_bsn_nfo!(bsn_nfo)
            return next_c
        end
        return 0
    elseif iseven(next_c) &&   (max_check <= bsn_nfo.consecutive_match < 2*max_check)
        # We make sure we hit the attractor 60 consecutive times
        bsn_nfo.consecutive_match+=1
        return 0
    elseif iseven(next_c) && bsn_nfo.consecutive_match >= max_check*2
        # We have checked the presence of an attractor: tidy up everything and get a new box.
        #ind = (bsn_nfo.basin .== bsn_nfo.current_color+1)
        #bsn_nfo.basin[ind] .= 1
        find_and_replace!(bsn_nfo, bsn_nfo.current_color+1, 1)

        #bsn_nfo.basin[n] = bsn_nfo.current_color
        # store color
        set_data!(bsn_nfo.basin, n, bsn_nfo.current_color)
        store_attractor!(bsn_nfo, u) # store attractor
        bsn_nfo.current_color = bsn_nfo.next_avail_color
        bsn_nfo.next_avail_color += 2
        println("even and > max check ", next_c)
        reset_bsn_nfo!(bsn_nfo)
        return next_c+1;
    else
        return 0
    end
end


function store_attractor!(bsn_nfo::BasinInfo, u)
    if haskey(bsn_nfo.attractors , bsn_nfo.current_color)
        push!(bsn_nfo.attractors[bsn_nfo.current_color],  u) # store attractor
    else
        bsn_nfo.attractors[bsn_nfo.current_color] = Dataset([SVector(u[1],u[2])])  # init dic
    end
end


"""
    draw_basin!(xg, yg, integ, iter_f!::Function, reinit_f!::Function)
Compute an estimate of the basin of attraction on a two-dimensional plane. This is a low level function,
for higher level functions see: `basin_poincare_map`, `basin_discrete_map`, `basin_stroboscopic_map`

## Arguments
* `xg`, `yg` : 1-dim range vector that defines the grid of the initial conditions to test.
* `integ` : integrator handle of the dynamical system.
* `iter_f!` : function that iterates the map or the system, see step! from DifferentialEquations.jl and
examples for a Poincaré map of a continuous system.
* `reinit_f!` : function that sets the initial condition to test on a two dimensional projection of the phase space.
"""
function draw_basin!(xg, yg, integ, iter_f!::Function, reinit_f!::Function, get_u::Function, Ncheck)

    complete = 0;

    bsn_nfo = BasinInfo(ones(Int16, length(xg), length(yg)), xg, yg, iter_f!, reinit_f!, get_u, 2,4,0,0,0,1,1,0,0,Dict{Int16,Dataset{2,Float64}}(),0,(),())

    reset_bsn_nfo!(bsn_nfo)

    I = CartesianIndices(bsn_nfo.basin)
    j  = 1

    while complete == 0
        # pick the first empty box
        ind = 0
        for k in j:length(bsn_nfo.basin)
            if bsn_nfo.basin[I[k]] == 1
                j = k
                ind=I[k]
                break
            end
        end

        if ind == 0
            # We are done
            complete = 1
            break
        end

         ni = ind[1]; mi = ind[2];
         x0 = xg[ni]; y0 = yg[mi]

         # Tentatively assign a color: odd is for basins, even for attractors.
         # First color is one
         bsn_nfo.basin[ni,mi] = bsn_nfo.current_color + 1

         u0=[x0, y0]

         bsn_nfo.basin[ni,mi] = get_color_point!(bsn_nfo, integ, u0; Ncheck=Ncheck)
    end

    bsn_nfo.Na = Int((bsn_nfo.current_color-2)/2)

    return bsn_nfo
end



function get_color_point!(bsn_nfo::BasinInfo, integ, u0; Ncheck=2)
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

       n = get_box(new_u, bsn_nfo)

       if !isnothing(n) # apply procedure only for boxes in the defined space
           done = procedure!(bsn_nfo, n, new_u, Ncheck)
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




function get_box(u, bsn_nfo::BasinInfo)
    xg = bsn_nfo.xg; yg = bsn_nfo.yg; # aliases
    xstep = (xg[2]-xg[1])
    ystep = (yg[2]-yg[1])
    xu=u[1]
    yu=u[2]
    if typeof(bsn_nfo.basin) <: Cell
        if xu >= xg[1] && xu <= xg[end] && yu >= yg[1] && yu <= yg[end]
            return findleaf(bsn_nfo.basin,u)
        else
            return nothing
        end
    end


    n=0; m=0;
    # check boundary
    if xu >= xg[1] && xu <= xg[end] && yu >= yg[1] && yu <= yg[end]
        n = Int(round((xu-xg[1])/xstep))+1
        m = Int(round((yu-yg[1])/ystep))+1 # +1 for 1 based indexing
    else
        return nothing
    end
    return CartesianIndex(n,m)
end


function check_outside_the_screen!(bsn_nfo::BasinInfo, new_u, old_u, inlimbo)

    if norm(new_u-old_u) < 1e-5
        #println("Got stuck somewhere, Maybe an attractor outside the screen: ", new_u)
        #ind = (bsn_nfo.basin .== bsn_nfo.current_color+1)
        #bsn_nfo.basin[ind] .= 1
        find_and_replace!(bsn_nfo,bsn_nfo.current_color+1, 1)
        reset_bsn_nfo!(bsn_nfo)
        # this CI goes to a attractor outside the screen, set to -1 (even color)
        return -1  # get next box
    elseif inlimbo > 60*20
        #println("trajectory diverges: ", new_u)
        #ind = (bsn_nfo.basin .== bsn_nfo.current_color+1)
        #bsn_nfo.basin[ind] .= 1
        find_and_replace!(bsn_nfo,bsn_nfo.current_color+1, 1)
        reset_bsn_nfo!(bsn_nfo)
        # this CI is problematic or diverges, set to -1 (even color)
        return -1  # get next box
    end
    return 0
end


function find_and_replace!(bsn_nfo::BasinInfo, old_c, new_c)
    if typeof(bsn_nfo.basin) <: Matrix
        ind = (basin .== old_c)
        basin[ind] .= new_c
    else
        # while !isempty(bsn_nfo.visited)
        #     point = pop!(bsn_nfo.visited)
        #     lf = findleaf(bsn_nfo.basin, point)
        #     if lf.data == old_c
        #         lf.data = new_c
        #         if new_c == 1
        #             push!(bsn_nfo.available,point)
        #         end
        #     end
        # end
        for lf in allleaves(bsn_nfo.basin)
            if lf.data==old_c
                lf.data==new_c
            end
        end

    end
end




function reset_bsn_nfo!(bsn_nfo::BasinInfo)
    #@show bsn_nfo.step
    bsn_nfo.consecutive_match = 0
    bsn_nfo.consecutive_other_basins = 0
    bsn_nfo.prevConsecutives = 0
    bsn_nfo.prev_attr = 1
    bsn_nfo.prev_bas = 1
    bsn_nfo.prev_step = 0
    bsn_nfo.step = 0
end

function basins_map2D_tree(xg, yg, integ; T=0., Ncheck = 2)
    if T>0
        iter_f! = (integ) -> step!(integ, T, true)
    else
        iter_f! = (integ) -> step!(integ)
    end
    reinit_f! =  (integ,y) -> reinit!(integ, y)
    get_u = (integ) -> integ.u

    return draw_basin_tree!(xg, yg, integ, iter_f!, reinit_f!, get_u, Ncheck)
end


import RegionTrees: AbstractRefinery, needs_refinement, refine_data

struct MyRefinery <: AbstractRefinery
    tolerance::Float64
    init_dat::Bool
end

# These two methods are all we need to implement
function needs_refinement(r::MyRefinery, cell)
    if r.init_dat
        return maximum(cell.boundary.widths) > r.tolerance
    else
        @show cell
        return maximum(cell.boundary.widths) > r.tolerance #&& (length(unique([l.data for l in children(cell.parent)])) > 1)
    # else
    #     return false
    end
end

function refine_data(r::MyRefinery, cell::Cell, indices)
        return 1
end

function is_in_cell(cell::Cell,  point::AbstractVector, indices)
    or = cell.boundary.origin
    wd = cell.boundary.widths
    dv = cell.divisions
    if point[1] > or[1] + wd[1] || point[1] < or[1]
        return false
    elseif point[2] > or[2] + wd[2] || point[2] < or[2]
        return false
    else
        ind = [1,1]
        if point[1] > dv[1]
            ind[1]=2
        end
        if point[2] > dv[2]
            ind[2]=2
        end
        if ind[1] == indices[1] && ind[2] == indices[2]
            return true
        else
            return false
        end
    end
end


function get_att_index(cell::Cell,  point::AbstractVector)
    or = cell.boundary.origin
    wd = cell.boundary.widths
    dv = cell.divisions
    ind = 1,1
    if point[1] > dv[1]
        ind[1]=2
    end

    if point[2] > dv[2]
        ind[2]=2
    end

    return ind
end

function refine_data(r::MyRefinery, cell::Cell, indices, bsn_nfo::BasinInfo)
    if iseven(cell.data)
        att = cell.data
        for p in bsn_nfo.attractors[att]
            if is_in_cell(cell,p, indices)
                @show att,p,cell.boundary,indices
                return att
            end
        end
    end
    push!(bsn_nfo.available,cell.divisions)
    return 1
end

function get_neighbor_cells(root::Cell,cell::Cell)
    or = cell.boundary.origin
    wd = cell.boundary.widths
    dv = cell.divisions
    W=findleaf(root,dv + SVector(wd[1],0)).data
    E=findleaf(root,dv .+ SVector(-wd[1],0)).data
    N=findleaf(root,dv .+ SVector(0,wd[2])).data
    S=findleaf(root,dv .+ SVector(0,-wd[2])).data
    #@show [W E N S cell.data]
    return length(unique([W E N S cell.data])) > 1
end

function check_if_complete!(root::Cell, refinery::AbstractRefinery, bsn_nfo::BasinInfo)
    complete = true
    refine_function = (cell, indices) -> refine_data(refinery, cell, indices, bsn_nfo)
    leaf_stack = Dataset([0. 0.])
    for lf in allleaves(root)
        #if maximum(lf.boundary.widths) > refinery.tolerance && (length(unique([l.data for l in children(lf.parent)])) > 1)
        if maximum(lf.boundary.widths) > refinery.tolerance && get_neighbor_cells(root,lf)
            complete = false
            #split!(lf, refine_function)
            push!(leaf_stack,lf.boundary.origin)
        end
    end

    for lf in leaf_stack[2:end]
        l = findleaf(root,lf)
        split!(l, refine_function)
    end

    return complete
end



function draw_basin_tree!(xg, yg, integ, iter_f!::Function, reinit_f!::Function, get_u::Function, Ncheck)

    complete = 0;

    xi=xg[1]; yi=yg[1]; xf=xg[end]; yf=yg[end]

    bsn_nfo = BasinInfo(Cell(SVector(xi, yi), SVector(xf-xi, yf-xi),1), xg, yg, iter_f!, reinit_f!, get_u, 2,4,0,0,0,1,1,0,0,Dict{Int16,Dataset{2,Float64}}(),0,Deque{Vector{Float64}}(),Deque{Vector{Float64}}())

    reset_bsn_nfo!(bsn_nfo)

    # Initial refinement of the grid.
    r = MyRefinery(0.2,true)
    @show bsn_nfo.basin
    adaptivesampling!(bsn_nfo.basin, r)

    for lf in allleaves(bsn_nfo.basin)
        push!(bsn_nfo.available,lf.divisions)
    end
    r_more = MyRefinery(0.1,false)

    while complete == 0
        # pick the first empty box
        # lf = nothing
        # for leaf in allleaves(bsn_nfo.basin)
        #     if leaf.data == 1
        #         lf=leaf
        #     end
        # end
        #
        # if isnothing(lf)
        lf = nothing
        while !isempty(bsn_nfo.available)
            point = pop!(bsn_nfo.available)
            lf = findleaf(bsn_nfo.basin,point)
            if lf.data == 1
                break
            end
        end

        if isempty(bsn_nfo.available)
            println("refinement!")
            if check_if_complete!(bsn_nfo.basin, r_more, bsn_nfo)
                complete=1
                break
            end
            # for leaf in allleaves(bsn_nfo.basin)
            #     if leaf.data == 1
            #         lf=leaf
            #     end
            # end
        end
        # Get next available candidate
        #lf = findleaf(bsn_nfo.basin, pop!(bsn_nfo.available))


        # Tentatively assign a color: odd is for basins, even for attractors.
        # First color is one
        lf.data = bsn_nfo.current_color + 1

        u0 = RegionTrees.center(lf)

        lf.data = get_color_point!(bsn_nfo, integ, u0; Ncheck=Ncheck)

    end

    bsn_nfo.Na = Int((bsn_nfo.current_color-2)/2)

    return bsn_nfo
end

function get_data(bsn::Matrix,n)
    return bsn[n]
end

function get_data(bsn::Cell,n)
    return n.data
end

function set_data!(bsn::Matrix,n,c)
    bsn[n]=c
end

function set_data!(bsn::Cell,n,c)
    n.data=c
end
