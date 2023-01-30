struct SpectralVars{FT,C} <: AbstractVars
     u :: FT
     v :: FT
     ϕ :: FT
     c :: FT
    uh :: C
    vh :: C
    ϕh :: C
    ch :: C
end

struct FiniteDiffVars{FT} <: AbstractVars
    u :: FT
    v :: FT
    ϕ :: FT
    c :: FT
end

function SpectralVars(G::OneDGrid)

    T = eltype(G)
    Dev = typeof(G.device)

    @devzeros Dev T G.nx u v ϕ c
    @devzeros Dev Complex{T} G.nkr uh vh ϕh ch

    return SpectralVars(u, v, ϕ, c, uh, vh, ϕh, ch)

end

function FiniteDiffVars(G::OneDGrid)

    T = eltype(G)

    return FiniteDiffVars(
        zeros(T,grid.nx),
        zeros(T,grid.nx),
        zeros(T,grid.nx),
        zeros(T,grid.nx),
    )

end

function SpectralVars(G::TwoDGrid)

    T = eltype(G)
    Dev = typeof(G.device)

    @devzeros Dev T (G.nx, G.ny) u v ϕ c
    @devzeros Dev Complex{T} (G.nkr, G.nl) uh vh ϕh ch

    return SpectralVars(u, v, ϕ, c, uh, vh, ϕh, ch)

end

function FiniteDiffVars(G::TwoDGrid)

    T = eltype(G)

    return FiniteDiffVars(
        zeros(T,G.nx,G.ny),
        zeros(T,G.nx,G.ny),
        zeros(T,G.nx,G.ny),
        zeros(T,G.nx,G.ny),
    )

end