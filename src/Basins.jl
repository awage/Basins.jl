module Basins

#using DifferentialEquations
using DynamicalSystems
using LinearAlgebra
using ImageFiltering
using NearestNeighbors
using Combinatorics
using Roots

export compute_saddle, detect_wada_grid_method, detect_wada_merge_method
include("wada_detection.jl")
include("straddle.jl")

end
