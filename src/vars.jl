struct SpectralVars{FT,C} <: AbstractVars
      u :: FT
      v :: FT
      ϕ :: FT
      c :: FT
     ϕf :: FT
     uh :: C
     vh :: C
     ϕh :: C
     ch :: C
    ϕfh :: C
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

    @devzeros Dev T G.nx u v ϕ c ϕf
    @devzeros Dev Complex{T} G.nkr uh vh ϕh ch ϕfh

    return SpectralVars(u, v, ϕ, c, ϕf, uh, vh, ϕh, ch, ϕfh)

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

    @devzeros Dev T (G.nx, G.ny) u v ϕ c ϕf
    @devzeros Dev Complex{T} (G.nkr, G.nl) uh vh ϕh ch ϕfh

    return SpectralVars(u, v, ϕ, c, ϕf, uh, vh, ϕh, ch, ϕfh)

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