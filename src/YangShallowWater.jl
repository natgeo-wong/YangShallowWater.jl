module YangShallowWater

export GenerateGrid

## Modules Used
using Dates
using LinearAlgebra: mul!, ldiv!
using NCDatasets
using Reexport
using Printf

import Base: eltype, show

@reexport using FourierFlows

## Including other files in the module

include("setup.jl")
include("model.jl")

end
