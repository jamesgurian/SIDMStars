
"""
    SidmStars



"""

module SIDMStars
using LinearAlgebra
using DifferentialEquations
using HDF5
using ForwardDiff
using DataInterpolations
using SpecialFunctions
using QuadGK
using FunctionWrappers
import FunctionWrappers:FunctionWrapper    
using Unitful
using UnitfulAstro
using StaticArrays
using PyPlot
include("bvp.jl")
include("dm_params.jl")
include("transport.jl")
include("util.jl")
include("init.jl")
include("snapshot.jl")
include("evolve.jl")

module Constants
        using PhysicalConstants
        using ..SIDMStars: @import_everything
            @import_everything PhysicalConstants.CODATA2018
    end

export RunParams
export initialize_profile
export Constants
export nfw_profile
export dk_profile
export setup_dk_profile
export evolve
export read_snapshot_file
# Include implementation files (place implementation in src/bvp.jl)


# Export the public API



end # module SidmStars