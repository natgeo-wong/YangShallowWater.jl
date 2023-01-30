struct Convection1DGrid{FT<:Real} <: AbstractGrid
    nx :: Int
    cx :: Vector{Int}
    Fc :: Vector{FT}
end

struct Convection1DFlux{FT<:Real} <: AbstractGrid
     old :: Vector{Int}
     new :: Vector{Int}
      Δt :: Vector{FT}
end

struct Convection2DGrid{FT<:Real} <: AbstractGrid
    nx :: Int
    ny :: Int
    cx :: Vector{Int,2}
    cy :: Vector{Int,2}
    Fc :: Array{FT,2}
end

struct Convection2DFlux{FT<:Real} <: AbstractGrid
     old :: Array{Int,2}
     new :: Array{Int,2}
      Δt :: Array{FT,2}
end

struct Convection1D{FT<:Real} <: Forcing1D
    τc :: FT # Convective Damping Timescale
    ϕc :: FT # Convection-triggering geopotential
    rc :: FT # Convection radius
    Sc :: FT # Number DENSITY of convective events
    Fc :: Vector{FT}
    ConvectionGrid :: Convection1DGrid{FT}
    ConvectionFlux :: Convection1DFlux{FT}
    Equation :: Function
end

struct Convection2D{FT<:Real} <: Forcing2D
    τc :: FT # Convective Damping Timescale
    ϕc :: FT # Convection-triggering geopotential
    rc :: FT # Convection radius
    Sc :: FT # Number DENSITY of convective events
    Fc :: Array{FT,2}
    ConvectionGrid :: Convection2DGrid{FT}
    ConvectionFlux :: Convection2DFlux{FT}
    Equation :: Function
end

function GenerateConvection(
	τc :: Real,
    ϕc :: Real,
    rc :: Real,
    Sc :: Real,
    convectionfunction :: Function,
    G :: OneDGrid,
)

    T = eltype(G)

    r2 = rc^2
    dx = G.Lx / G.nx; nx_c = floor(rc/dx); xc = nx_c+1

    cflux = ones(T, 2*nx_c+1) ./ P.τc / P.rc
    for ix = -nx_c : nx_c
        cflux[xc+ix] *= (1 - (dx*ix)^2/r2)
    end

    return Convection1D{T}(
        τc, ϕc, rc, Sc,
        Convection1DGrid{T}(
            2*nx_c+1,
            collect(-nx_c:nx_c)
            cflux
        ),
        Convection1DFlux{T}(
            zeros(Int, G.nx),
            zeros(Int, G.nx),
            zeros(T, G.nx)
        ),
        convectionfunction
    )

end

function GenerateConvection(
	τc :: Real,
    ϕc :: Real,
    rc :: Real,
    Sc :: Real,
    convectionfunction :: Function,
    G :: TwoDGrid,
)

    T = eltype(G)

    r2 = rc^2
    dx = G.Lx / G.nx; nx_c = floor(rc/dx); xc = nx_c+1
    dy = G.Ly / G.ny; ny_c = floor(rc/dy); yc = ny_c+1

    cflux = ones(T, 2*nx_c+1, 2*ny_c+1) ./ P.τc / pi*rc2
    for iy = -ny_c : ny_c, ix = -nx_c : nx_c
        ir = sqrt((dx*ix)^2+(dy*iy)^2)
        cflux[xc+ix,yc+iy] *= (1 - ir/r2)
    end

    return Convection2D{T}(
        τc, ϕc, rc, Sc,
        Convection2DGrid{T}(
            2*nx_c+1, 2*ny_c+1,
            collect(-nx_c:nx_c),
            collect(-ny_c:ny_c),
            cflux
        ),
        Convection2DFlux{T}(
            zeros(Int, G.nx, G.ny),
            zeros(Int, G.nx, G.ny),
            zeros(T, G.nx, G.ny)
        ),
        convectionfunction
    )

end