module Basins

using DifferentialEquations
using LinearAlgebra
using LsqFit

export draw_basin, basin_entropy, uncertainty_dimension
# Write your package code here.
include("compute_ba.jl")
include("compute_be.jl")
include("compute_uncertain_dim.jl")

end
