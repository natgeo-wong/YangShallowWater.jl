struct ϕForcing1D{FT<:Real} <: Forcing1D
     τl :: FT   # Large-scale relaxation timescale
     ϕ0 :: FT   # Large-scale geopotential
     ϕl :: FT   # Large-scale geopotential forcing
    wtg :: Bool # Is the Weak-Temperature Gradient scheme used?
    ϕForcingGrid :: Vector{FT}
    Equation :: Function
end

struct ϕForcing2D{FT<:Real} <: Forcing2D
     τl :: FT   # Large-scale relaxation timescale
     ϕ0 :: FT   # Large-scale geopotential
     ϕl :: FT   # Large-scale geopotential forcing
    wtg :: Bool # Is the Weak-Temperature Gradient scheme used?
    ϕForcingGrid :: Array{FT,2}
    Equation :: Function
end

function GenerateϕForcing(
	τl  :: Real,
    ϕ0  :: Real,
    ϕl  :: Real,
    wtg :: Bool,
    ϕfunction :: Function,
    G :: OneDGrid,
)

    T = eltype(G)

    ϕflux = ones(T, G.nx)

    return ϕForcing1D{T}(τl, ϕ0, ϕl, wtg, ϕflux, ϕfunction)

end

function GenerateϕForcing(
	τl  :: Real,
    ϕ0  :: Real,
    ϕl  :: Real,
    wtg :: Bool,
    ϕfunction :: Function,
    G :: TwoDGrid,
)

    T = eltype(G)

    ϕflux = ones(T, G.nx, G.ny)

    return ϕForcing2D{T}(τl, ϕ0, ϕl, wtg, ϕflux, ϕfunction)

end