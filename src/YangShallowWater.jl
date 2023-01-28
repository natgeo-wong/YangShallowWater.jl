module YangShallowWater

## Modules Used
using Dates
using FourierFlows
using LinearAlgebra: mul!, ldiv!
using NCDatasets
using Printf
import Base: eltype, show



## Including other files in the module

include("setup.jl")
include("model.jl")

end
