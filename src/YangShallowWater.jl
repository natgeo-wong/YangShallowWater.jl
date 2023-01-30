module YangShallowWater

## Modules Used
using Dates
using LinearAlgebra: mul!, ldiv!
using NCDatasets
using Reexport
using Printf

import Base: eltype, show, run

@reexport using FourierFlows

export
        GenerateGrid, DefineParams, CreateModel,
        run

## Abstract SuperTypes
"""
    AbstractModel

Abstract supertype for different models that YangShallowWater can run.
"""
abstract type AbstractModel end

"""
    YSWParams

Abstract supertype for different parameter inputs (i.e., forcing, convection, etc.)
"""
abstract type YSWParams <: AbstractParams end
abstract type ForcingParams <: YSWParams end

"""
    AbstractForcing

Abstract supertype for different parameter inputs (i.e., forcing, convection, etc.)
"""
abstract type AbstractForcing end

"""
    YSWParams

Abstract supertype for different parameter inputs (i.e., forcing, convection, etc.)
"""
abstract type Forcing1D <: AbstractForcing end
abstract type Forcing2D <: AbstractForcing end

## Including other files in the module

include("grid.jl")

include("param/convection.jl")
include("param/potentialforcing.jl")
include("param/define.jl")

include("vars.jl")

include("model.jl")

include("calculate/convection.jl")
include("calculate/potential.jl")

include("equations/spectral2D.jl")

include("run.jl")

end
