
function get_boundary_filt(basin)

    # Kernel
    w = centered([1 1 1; 1 -8 1 ; 1 1 1])
    # replicate for boundary conditions
    res =imfilter(basin,w,"replicate")
    # segmentation
    return abs.(res) .> 0

end

function get_list(basin,y1 ,y2)

    num_att = length(unique(basin))

    # v_list empty container that will receive the different coordinates corresponding to each basins,
    # Each element of this set is an array of coordinates of different sizes.
    v_list = []

    # original boundary
    bnd = get_boundary_filt(basin)
    I1=findall(bnd .== 1);
    push!(v_list , hcat([[y1[ind[1]]; y2[ind[2]]] for ind in I1 ]...))

    # Merging !!
    basin_i_j = zeros(Int8,size(basin))
    for k in unique(basin)
        fill!(basin_i_j, zero(0)) # erase matrix
        I = findall(basin .== k)
        [basin_i_j[x] = 1 for x in I];
        # compute boundary
        bnd = get_boundary_filt(basin_i_j)
        I1=findall(bnd .== 1);
        # Get coordinates lists from matrix
        push!(v_list , hcat([[y1[ind[1]]; y2[ind[2]]] for ind in I1 ]...))
    end
    return v_list
end

function wada_merge_dist(basin,y1,y2)
    num_att = length(unique(basin))

    v_list=get_list(basin, y1, y2)

    # compute distances using combbinations of 2 elements from a collection
    ind = combinations(1:num_att,2)

    min_dist = Inf;
    max_dist = 0
    for (k, patron) in enumerate(ind)

       v01 = v_list[patron[1]]
       v02 = v_list[patron[2]]

       hd = haussdorff_dist(v01,v02)

       if hd  > max_dist
           max_dist = hd
       end

       if hd  < min_dist
           min_dist = hd
       end

  end

   return max_dist, min_dist
end

function haussdorff_dist(v01,v02)
       kdtree1 = KDTree(v01)
       idxs12, dists12 = knn(kdtree1, v02, 1, true)
       kdtree2 = KDTree(v02)
       idxs21, dists21 = knn(kdtree2, v01, 1, true)
       max12 = maximum(dists12);
       max21 = maximum(dists21);
       hd =max(max12[1],max21[1])
       return hd
end




mutable struct ode_info{I}
    bsn_nfo::basin_info # basin info for BA routine
    integ::I               # integrator
end


function reset_point_data!(ode_nfo)
    reset_bsn_nfo!(ode_nfo.bsn_nfo)
end


function detect_wada_grid_method(integ, bsn_nfo; max_iter=10)

   ode_nfo = ode_info(bsn_nfo, integ)
   num_att = length(unique(bsn_nfo.basin))

   # helper function to obtain coordinates
   index_to_coord(p) = [bsn_nfo.xg[p[1]], bsn_nfo.yg[p[2]]]

   # obtain points in the boundary
   bnd = get_boundary_filt(bsn_nfo.basin)
   p1_ind = findall(bnd .> 0)

   # initialize empty array of indices and collection of empty sets of colors
   clr_mat = [Set{Int16}() for i=1:length(p1_ind)]
   p2_ind = typeof(p1_ind)(undef, length(p1_ind))
   W = zeros(num_att, max_iter)

   # Initialize matrices (step 1)
   for (k,p1) in enumerate(p1_ind)
       p2, nbgs = get_neighbor_and_colors(bsn_nfo.basin, [p1[1], p1[2]])
       if length(nbgs) > 1
           # keep track of different colors and neighbor point
           push!(clr_mat[k],nbgs...)
           p2_ind[k] = p2
       else
           println("Not in the boundary: weird...")
       end
       W[length(clr_mat[k]),1] += 1
   end

   # Do the iteration!
   for n = 2:max_iter
       for k in 1:length(p1_ind)
           pc1 = index_to_coord(p1_ind[k])
           pc2 = index_to_coord(p2_ind[k])
           # update number of colors
           clr_mat[k]=divide_and_test_W(ode_nfo, pc1, pc2, n, clr_mat[k], num_att)
           # update W matrix
           W[length(clr_mat[k]),n] += 1
       end
       @show W

       # Stopping criterion: if W[Na] in % increases less than ε  then stop or if we have more than 95% boxes in Na
       if  (abs(W[num_att,n] - W[num_att,n-1]) <0.01*W[num_att,n] || W[num_att,n] == 0.0 || W[num_att,n]/sum(W[:,n])>0.95) &&  n > num_att # make at least num_att iterations for stats
           return W[:,n]./sum(W[:,1])
       end

   end

   return W[:,end]./sum(W[:,1])
end


function get_neighbor_and_colors(basin, p)
    radius = 1
    n=p[1]
    m=p[2]
    v=Int16[]
    p2=CartesianIndex(-1,-1)
    # check neihbors and collect basin colors
    for k=n-radius:n+radius, l=m-radius:m+radius
        try
            push!(v,basin[k,l])
            if k != n || l != m
                if basin[n,m] != basin[k,l]
                    p2=CartesianIndex(k,l)
                end
            end
        catch
            #println("error in push!")
        end
    end
    return p2,unique(v)
end



function divide_and_test_W(ode_nfo, p1, p2, nstep, clrs, Na)

    if length(clrs)  == Na
        return clrs
    end

    lerp(t) = (1-t)*p1 + t*p2 # linear interpolation function

    # generates coordinates for initial conditions to test with linear interpolation
    npts = 2^(nstep)
    vecs = Set([lerp(ratio) for ratio in range(0,1, step=1/npts)])  # create a set with the points interspeded
    npts = 2^(nstep-1)
    vecs_before = Set([lerp(ratio) for ratio in range(0,1, step=1/npts)])  # remove points already computed in the previous step
    pnt_to_test = setdiff(vecs,vecs_before) # get the points we have to test

    # get colors and update color set for this box!
    for pnt in pnt_to_test
        clr = get_color_point!(ode_nfo.bsn_nfo, ode_nfo.integ, pnt[1],pnt[2])
        push!(clrs, clr)
        if length(clrs)  == Na
            break
        end
    end

    return clrs

end
