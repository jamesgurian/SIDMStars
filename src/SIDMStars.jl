
"""
    SidmStars

A lightweight Julia package layout for SIDM star structure and boundary-value problems.


"""

module SIDMStars
using LinearAlgebra
using DifferentialEquations
import BoundaryValueDiffEq as BVP
using Unitful
using UnitfulAstro
using SpecialFunctions
using QuadGK
using FunctionWrappers
import FunctionWrappers:FunctionWrapper    


include("bvp.jl")
include("dm_params.jl")
include("transport.jl")
include("util.jl")
include("init.jl")
module Constants
        using PhysicalConstants
        using ..SIDMStars: @import_everything
            @import_everything PhysicalConstants.CODATA2018
    end

export RunParams
export initialize_profile
export Constants

# Include implementation files (place implementation in src/bvp.jl)


# Export the public API



end # module SidmStars