struct Model1DSpectral{FT} <: AbstractModel
    Grid       :: OneDGrid
    Parameters :: YSWParams
    Variables  :: SpectralVars
end

struct Model1DFiniteDiff{FT} <: AbstractModel
    Grid       :: OneDGrid
    Parameters :: YSWParams
    Variables  :: FiniteDiffVars
end

struct Model2DSpectral{FT} <: AbstractModel
    Grid       :: TwoDGrid
    Parameters :: YSWParams
    Variables  :: SpectralVars
end

struct Model2DFiniteDiff{FT} <: AbstractModel
    Grid       :: TwoDGrid
    Parameters :: YSWParams
    Variables  :: FiniteDiffVars
end

function CreateModel(
    G :: OneDGrid,
    P :: YSWParams;
    spectral :: Bool = true,
)

    T = eltype(G)

    if spectral
          return Model1DSpectral{T}(G,P,SpectralVars())
    else; return Model1DFiniteDiff{T}(G,P,FiniteDiffVars())
    end

end

function CreateModel(
    G :: TwoDGrid,
    P :: YSWParams;
    spectral :: Bool = true,
)

    T = eltype(G)

    if spectral
          return Model2DSpectral{T}(G,P,SpectralVars())
    else; return Model2DFiniteDiff{T}(G,P,FiniteDiffVars())
    end

end