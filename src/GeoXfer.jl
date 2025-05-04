module GeoXfer

using LinearAlgebra
using GLMakie
using JuMP, Ipopt
using StaticArrays
using DifferentialEquations
using Optim
using FastGaussQuadrature
using JLD2, FileIO
using LegendrePolynomials

include("Utils.jl")
include("Dynamics.jl")
include("TypeDefs.jl")
include("UnivKepler.jl")
include("BiHohmann.jl")
include("BasicLambert.jl")
include("Primers.jl")
include("Plotting.jl")
include("Pseudospectral.jl")
include("PseudospectralImpulsive.jl")

export basic_lambert

const DISTANCE_UNIT = 6378.1363
const TIME_UNIT = 806.81099130673
export DISTANCE_UNIT, TIME_UNIT

end # module GeoXfer
