module Basins

#using DifferentialEquations
using DynamicalSystems
using LinearAlgebra
using ImageFiltering
using NearestNeighbors
using Combinatorics
using Roots
using LsqFit

export basin_entropy, detect_wada_grid_method, detect_wada_merge_method
export box_counting_dim
export draw_basin!, basin_map, basin_general_ds
export basin_stability, uncertainty_exponent

export compute_saddle,compute_basin_precise
# Write your package code here.
include("compute_ba.jl")
include("compute_be.jl")
include("compute_uncertain_dim.jl")
include("compute_wada.jl")
include("straddle.jl")

end
