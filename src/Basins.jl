module Basins

#using DifferentialEquations
using DynamicalSystems
using LinearAlgebra
using ImageFiltering
using NearestNeighbors
using Combinatorics
using Roots
using LsqFit
using RegionTrees
using DataStructures

export basin_entropy, detect_wada_grid_method, detect_wada_merge_method
export box_counting_dim
export draw_basin!, basin_map, basin_general_ds
export basin_stability, uncertainty_exponent, static_estimate

export compute_saddle,compute_basin_precise
# Write your package code here.
include("basins_attraction.jl")
include("basin_entropy.jl")
include("uncertain_dim.jl")
include("wada_detection.jl")
include("straddle.jl")

end
