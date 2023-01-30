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
    :: OneDGrid,
)

    T = eltype(G)

    ϕflux = ones(T, 2*nx_c+1) ./ P.τc / P.rc
    for ix = -nx_c : nx_c
        ϕflux[xc+ix] *= (1 - (dx*ix)^2/r2)
    end

    return ϕForcing1D{T}(τl, ϕ0, ϕl, wtg, ϕflux, ϕfunction)

end

function GenerateϕForcing(
	τl  :: Real,
    ϕ0  :: Real,
    ϕl  :: Real,
    wtg :: Bool,
    ϕfunction :: Function,
    :: TwoDGrid,
)

    T = eltype(G)

    ϕflux = ones(T, 2*nx_c+1, 2*ny_c+1) ./ P.τc / pi*rc2
    for iy = -ny_c : ny_c, ix = -nx_c : nx_c
        ir = sqrt((dx*ix)^2+(dy*iy)^2)
        ϕflux[xc+ix,yc+iy] *= (1 - ir/r2)
    end

    return ϕForcing2D{T}(τl, ϕ0, ϕl, wtg, ϕflux, ϕfunction)

end

function calculateConvection()



end