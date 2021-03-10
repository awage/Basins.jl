module Basins

using DifferentialEquations
using LinearAlgebra


export draw_basin, basin_entropy
# Write your package code here.
include("compute_ba.jl")
include("compute_be.jl")


end
