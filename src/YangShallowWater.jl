module YangShallowWater

export GenerateGrid

## Modules Used
using Dates
using LinearAlgebra: mul!, ldiv!
using NCDatasets
using Reexport
using Printf

import Base: eltype, show, run

@reexport using FourierFlows

export
        YSWParams,
        GenerateGrid, DefineParams,
        run

## Including other files in the module

include("setup.jl")
include("model.jl")
include("run.jl")

end
