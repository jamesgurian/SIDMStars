module SIDMStars
using LinearAlgebra
using LoopVectorization
using DifferentialEquations
# DifferentialEquations v8 no longer re-exports the solver algorithms;
# Rodas5 (bvp.jl) now lives in the OrdinaryDiffEqRosenbrock sublibrary.
using OrdinaryDiffEqRosenbrock
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
using Roots
using PrecompileTools
using NumericalIntegration
using NonlinearSolve

using RecursiveArrayTools
include("common.jl")

include("bvp.jl")
include("transport.jl")
include("util.jl")
include("init.jl")
include("snapshot.jl")
include("convect.jl")
include("evolve.jl")



module Constants
        using PhysicalConstants
        using ..SIDMStars: @import_everything
            @import_everything PhysicalConstants.CODATA2018
    end

export RunParams
export CrossSection
export PowerLawCS
export InterpolatedCS
export RegularizedCS
export ConvectionOn
export ConvectionOff
export EvaporationOff
export EvaporationOn
export initialize_profile
export Constants
export nfw_profile
export dk_profile
export setup_dk_profile
export evolve
export read_snapshot_file


@setup_workload begin
    
    @compile_workload begin
        redirect_stdout(devnull) do
        include("precompile_load.jl")
        end
    end
end

end 