struct YSWParams{FT<:Real} <: AbstractParams
    c  :: FT # Gravity Wave speed
    τd :: FT # Linear Damping Timescale
    τc :: FT # Convective Damping Timescale
    τl :: FT # WTG relaxation timescale
    ϕ0 :: FT # Large-scale geopotential
    ϕc :: FT # Convection-triggering geopotential
    rc :: FT # Convection radius
    nc :: FT # Number DENSITY of convective events
end

function GenerateGrid(
	FT = Float64;
	nx :: Int,
    ny :: Int = 1,
    Lx :: FT,
    Ly :: FT = 0
)

    if isone(ny) || iszero(Ly)
        return OneDGrid(;
            nx=nx, Lx=Lx, x0=Lx/(2nx)
        )
    else
        return TwoDGrid(;
            nx, Lx, x0=Lx/(2nx),
            ny, Ly, y0=Ly/(2ny)
        )
    end

end

DefineParams(
	FT = Float64;
	c  :: FT = 20.,     # units in m s**-1
    τd :: FT = 0.6,     # units in hours
    τc :: FT = 1.,      # units in days
    τl :: FT = 1.,      # units in days
    ϕ0 :: FT = 0,       # units in m**2 s**-2
    ϕc :: FT = c^2,     # units in m**2 s**-2
    rc :: FT = 10.,     # units in km
    nc :: FT = 4.e-10,  # units in m**-2 s**-1
) = YSWParams{FT}(c, τd*3600, τc*86400, τl*86400, ϕ0, ϕc, rc*1000, nc)