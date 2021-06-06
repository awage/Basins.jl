"""
    basin_entropy(basin; eps_x=20, eps_y=20)
This algorithm computes the basin entropy of a computed basin of attraction on a regular grid.
The function return the basin entropy and the boundary basin entropy.

[A. Daza, A. Wagemakers, B. Georgeot, D. Guéry-Odelin and M. A. F. Sanjuán, Basin entropy:
a new tool to analyze uncertainty in dynamical systems, Sci. Rep., 6, 31416, (2016).]

## Arguments
* `basin` : the matrix containing the information of the basin.

## Keyword arguments
* `eps_x`, `eps_y` : define the size of the box that will be used to sample the basin
"""
function basin_entropy(basin; eps_x=20, eps_y=20)
    r,c= size(basin)
    vals = unique(basin)
    S=Int16(length(vals))
    pn=zeros(Float64,1,S)
    Sb=0
    Nb=0
    N=0
    for x in 1:eps_x:(r-eps_x+1), y in 1:eps_y:(c-eps_y+1)
        I = CartesianIndices((x:x+eps_x-1,y:y+eps_y-1))
        box_values=[basin[k] for k in I]
        N=N+1
        for (k,v) in enumerate(vals)
            pn[k]=count(x->x==v,box_values)/length(box_values)
        end
        Nb = Nb + (length(unique(box_values))>1)
        Sb = Sb + sum(entropy.(pn))
    end
        return Sb/N, Sb/Nb
end


function entropy(p::Float64)
    if p == 0.
        h = 0.
    else
        h=p*log(1/p)
    end
    return h
end

function basin_entropy(bsn::BasinInfo; eps_x=20, eps_y=20)
    ind  = findall(iseven.(bsn.basin) .== true)
    # Not sure it is necessary to deepcopy
    basin_test = deepcopy(bsn.basin)
    # Remove attractors from the picture
    for k in ind; basin_test[k] = basin_test[k]+1; end
    basin_entropy(basin_test; eps_x=eps_x, eps_y=eps_y)
end

"""
    basin_stability(basin)
This algorithm computes the basin stability of a computed basin of attraction on a regular grid.
This function returns a vector with the relative volume of the basins.

[P. Menck, J. Heitzig, N. Marwan et al. How basin stability complements the linear-stability paradigm. Nature Phys 9, 89–92 (2013).]

## Arguments
* `basin` : the matrix containing the information of the basin.

"""
function basin_stability(basin)

N = length(unique(basin))
v = zeros(1,N)

for (k,b) in enumerate(unique(basin))
        v[k] = count(basin .== b)/length(basin)
end

return v

end


function basin_stability(bsn::BasinInfo)
    ind  = findall(iseven.(bsn.basin) .== true)
    # Not sure it is necessary to deepcopy
    basin_test = deepcopy(bsn.basin)
    # Remove attractors from the picture
    for k in ind; basin_test[k] = basin_test[k]+1; end
    return basin_stability(basin_test)
end



"""
    fractal_test(bsn::BasinInfo)

## Arguments
* `bsn` :

## Keyword arguments
*

"""
function fractal_test(bsn::BasinInfo)

# disk radius
εₕ=5

# Sanity check.
if length(bsn.xg)/εₕ < 50
    @warn "Maybe the size of the grid is not fine enough."
end

nx = length(bsn.xg)
ny = length(bsn.yg)
xi=bsn.xg[1]; xf=bsn.xg[end]; yi=bsn.yg[1]; yf=bsn.yg[end];


# Now get the values in the boxes.
while Nb < 1e4

x = rand()*(xf-xi)+xi
y = rand()*(yf-yi)+yi
ind = find_in_range(x,y,εₕ)

end


end

function find_in_range(x,y,εₕ)


end
