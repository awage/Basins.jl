module Basins

#using DifferentialEquations
using DynamicalSystems
using LinearAlgebra
using LsqFit
using ImageFiltering
using NearestNeighbors
using Combinatorics
using Roots

export draw_basin, basin_entropy, uncertainty_dimension_sample,compute_wada_W, wada_merge_dist
export basin_poincare_map, basin_stroboscopic_map, basin_discrete_map, box_counting_dim
# Write your package code here.
include("compute_ba.jl")
include("compute_be.jl")
include("compute_uncertain_dim.jl")
include("compute_wada.jl")

end
