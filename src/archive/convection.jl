struct Convection1D{FT<:Real} <: Forcing1D
    τc :: FT # Convective Damping Timescale
    ϕ0 :: FT # Large-scale geopotential
    ϕc :: FT # Convection-triggering geopotential
    rc :: FT # Convection radius
    Sc :: FT # Number DENSITY of convective events
    ConvectionFlux  :: Convection1DFlux{FT}
    ConvectionGrid  :: Convection1DGrid{FT}
    calcConvection! :: Function
end

struct Convection2D{FT<:Real} <: Forcing2D
    τc :: FT # Convective Damping Timescale
    ϕ0 :: FT # Large-scale geopotential
    ϕc :: FT # Convection-triggering geopotential
    rc :: FT # Convection radius
    Sc :: FT # Number DENSITY of convective events
    ConvectionFlux  :: Convection2DFlux{FT}
    ConvectionGrid  :: Convection2DGrid{FT}
    calcConvection! :: Function
end

struct Convection1DFlux{FT<:Real} <: AbstractGrid
    Fc :: Vector{FT}
end

struct Convection2DFlux{FT<:Real} <: AbstractGrid
    Fc :: Array{FT,2}
end

struct Convection1DGrid{FT<:Real} <: AbstractGrid
     old :: Vector{Int}
     new :: Vector{Int}
      Δt :: Vector{FT}
    grid :: Vector{FT}
end

struct Convection2DGrid{FT<:Real} <: AbstractGrid
     old :: Array{Int,2}
     new :: Array{Int,2}
      Δt :: Array{FT,2}
    grid :: Array{FT,2}
end

function GenerateConvection(
	P :: YSWParams{FT},
    G :: OneDGrid,
    F :: Function,
) where FT <: Real

    rc = P.rc
    r2 = P.rc^2
    dx = G.Lx / G.nx; nx_c = floor(rc/dx); xc = nx_c+1

    cflux = ones(FT, 2*nx_c+1) ./ P.τc / P.rc
    for ix = -nx_c : nx_c
        cflux[xc+ix] *= (1 - (dx*ix)^2/r2)
    end

    return Convection1D(
        P.τc, P.ϕ0, P.ϕc, P.rc, P.Sc,
        Convection1DFlux(cflux),
        Convection1DGrid(
            zeros(Int, G.nx),
            zeros(Int, G.nx),
            zeros(FT, G.nx),
            zeros(FT, G.nx)
        ),
        F
    )

end

function GenerateConvection(
	P :: YSWParams{FT},
    G :: TwoDGrid,
    F :: Function,
) where FT <: Real

    rc = P.rc
    r2 = P.rc^2
    dx = G.Lx / G.nx; nx_c = floor(rc/dx); xc = nx_c+1
    dy = G.Ly / G.ny; ny_c = floor(rc/dy); yc = ny_c+1

    cflux = ones(FT, 2*nx_c+1, 2*ny_c+1) ./ P.τc / pi*rc2
    for iy = -ny_c : ny_c, ix = -nx_c : nx_c
        ir = sqrt((dx*ix)^2+(dy*iy)^2)
        cflux[xc+ix,yc+iy] *= (1 - ir/r2)
    end

    return Convection1D(
        P.τc, P.ϕ0, P.ϕc, P.rc, P.Sc,
        Convection2DFlux(cflux),
        Convection2DGrid(
            zeros(Int, G.nx, G.ny),
            zeros(Int, G.nx, G.ny),
            zeros(FT, G.nx, G.ny),
            zeros(FT, G.nx, G.ny)
        ),
        F
    )

end

function calculateConvection()



end