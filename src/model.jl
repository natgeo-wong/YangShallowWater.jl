struct Model1DSpectral{FT} <: AbstractModel
    Grid       :: OneDGrid
    Parameters :: YSWParams
    Variables  :: SpectralVars
    Directory  :: String
end

struct Model1DFiniteDiff{FT} <: AbstractModel
    Grid       :: OneDGrid
    Parameters :: YSWParams
    Variables  :: FiniteDiffVars
    Directory  :: String
end

struct Model2DSpectral{FT} <: AbstractModel
    Grid       :: TwoDGrid
    Parameters :: YSWParams
    Variables  :: SpectralVars
    Directory  :: String
end

struct Model2DFiniteDiff{FT} <: AbstractModel
    Grid       :: TwoDGrid
    Parameters :: YSWParams
    Variables  :: FiniteDiffVars
    Directory  :: String
end

function CreateModel(
    G :: OneDGrid,
    P :: YSWParams;
    spectral :: Bool = true,
    path :: String = homedir()
)

    T = eltype(G)

    if spectral
          return Model1DSpectral{T}(G,P,SpectralVars(G),path)
    else; return Model1DFiniteDiff{T}(G,P,FiniteDiffVars(G),path)
    end

end

function CreateModel(
    G :: TwoDGrid,
    P :: YSWParams;
    spectral :: Bool = true,
    path :: String = homedir()
)

    T = eltype(G)

    if spectral
          return Model2DSpectral{T}(G,P,SpectralVars(G),path)
    else; return Model2DFiniteDiff{T}(G,P,FiniteDiffVars(G),path)
    end

end