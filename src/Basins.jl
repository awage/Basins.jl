module Basins

#using DifferentialEquations
using DynamicalSystems
using LinearAlgebra
using ImageFiltering
using NearestNeighbors
using Combinatorics
using Roots

export compute_saddle
include("wada_detection.jl")
include("straddle.jl")

end
