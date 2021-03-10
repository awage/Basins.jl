
#using StatsBase

# Compute Sb and Sbb Boundary Basin Entropy
function basin_entropy(cuenca, eps_x, eps_y)
r,c= size(cuenca)
vals = unique(cuenca)
S=Int16(length(vals))
pn=zeros(1,S)
Sb=0
Nb=0
N=0
for x = 1:eps_x:(r-eps_x+1)
    for y = 1:eps_y:(c-eps_y+1)
        x_coor=x:x+eps_x-1
        y_coor=y:y+eps_y-1
        box_values=[cuenca[k,m] for k in x_coor, m in y_coor]
        N=N+1
        for (k,v) in enumerate(vals)
            pn[k]=count(x->x==Float64(v),box_values)/length(box_values)
        end
        #push!(p0,pn[1])
        Nb = Nb + (length(unique(box_values))>1)
        Sb = Sb + sum(entropy.(pn))
    end
end
    return Sb/N, Sb/Nb
end

function entropy(p)
    if p == 0
        h = 0
    else
        h=p*log(1/p)
    end
    return h
end
