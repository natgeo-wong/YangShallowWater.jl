abstract type ShallowWaterData end

struct Variable1DData{FT<:Real}
    statistics :: Vector{FT}
    data :: Array{FT,2}
end

struct Variable2DData{FT<:Real}
    statistics :: Array{FT,2}
    data :: Array{FT,3}
end

struct ShallowWaterData1D{FT<:Real} <: ShallowWaterData
      t :: Vector{FT}
      ϕ :: Variable1DData{FT}
     ϕf :: Variable1DData{FT}
      u :: Variable1DData{FT}
      v :: Variable1DData{FT}
      c :: Variable1DData{FT}
     Fc :: Variable1DData{FT}
    cΔt :: Variable1DData{FT}
end

struct ShallowWaterData2D{FT<:Real} <: ShallowWaterData
      t :: Vector{FT}
      ϕ :: Variable2DData{FT}
     ϕf :: Variable2DData{FT}
      u :: Variable2DData{FT}
      v :: Variable2DData{FT}
      c :: Variable2DData{FT}
     Fc :: Variable2DData{FT}
    cΔt :: Variable2DData{FT}
end

function createsavedata(
    dt     :: Real,
    nsteps :: Int,
    nstats :: Int,
    nsave  :: Int,
    G :: OneDGrid
)

    T = eltype(G)

    if !iszero(mod(nsave,nstats))
        error("$(modulelog()) - nstats=$nstats is not a factor of nsave=$nsave")
    end

    totaltime = dt * nsteps
    statstime = dt * nstats
    savetime  = dt * nsave

    ndaysave  = Int(86400 / savetime)
    ntotal    = Int(totaltime/savetime)

    if ntotal > ndaysave
        return ndaysave, ShallowWater1DData{T}(
            zeros(T,ndaysave),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ndaysave)),
        )
    else
        return ntotal, ShallowWaterData1D{T}(
            zeros(T,ndaysave),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
            Variable1DData{T}(zeros(T,G.nx), zeros(T,G.nx,ntotal)),
        )
    end

end

function createsavedata(
    dt     :: Real,
    nsteps :: Int,
    nstats :: Int,
    nsave  :: Int,
    G :: TwoDGrid
)

    T = eltype(G)

    if !iszero(mod(nsave,nstats))
        error("$(modulelog()) - nstats=$nstats is not a factor of nsave=$nsave")
    end

    totaltime = dt * nsteps
    savetime  = dt * nsave

    ndaysave  = Int(86400 / savetime)
    ntotal    = Int(totaltime/savetime)

    if ntotal > ndaysave
        return ndaysave, ShallowWaterData2D{T}(
            zeros(T,ndaysave),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ndaysave)),
        )
    else
        return ntotal, ShallowWaterData2D{T}(
            zeros(T,ntotal),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
            Variable2DData{T}(zeros(T,G.nx,G.ny), zeros(T,G.nx,G.ny,ntotal)),
        )
    end

end

function updatestats!(
    data :: ShallowWaterData,
    prob,
    :: SimpleParams
)

    @. data.ϕ.statistics   += prob.vars.ϕ
    @. data.ϕf.statistics  += prob.vars.ϕf
    @. data.u.statistics   += prob.vars.u
    @. data.v.statistics   += prob.vars.v
    @. data.Fc.statistics  += prob.vars.c

end

function updatestats!(
    data :: ShallowWaterData,
    prob,
    :: ForcingParams
)

    @. data.ϕ.statistics   += prob.vars.ϕ
    @. data.ϕf.statistics  += prob.vars.ϕf
    @. data.u.statistics   += prob.vars.u
    @. data.v.statistics   += prob.vars.v
    @. data.c.statistics   += prob.params.convection.ConvectionFlux.new
    @. data.Fc.statistics  += prob.vars.c
    @. data.cΔt.statistics += prob.params.convection.ConvectionFlux.Δt

end

function savestats!(
    data :: ShallowWaterData1D,
    istats :: Int,
    itime :: Int,
)

    @views @. data.ϕ.data[:,itime]   = data.ϕ.statistics   / istats
    @views @. data.ϕf.data[:,itime]  = data.ϕf.statistics  / istats
    @views @. data.u.data[:,itime]   = data.u.statistics   / istats
    @views @. data.v.data[:,itime]   = data.v.statistics   / istats
    @views @. data.c.data[:,itime]   = data.c.statistics   / istats
    @views @. data.Fc.data[:,itime]  = data.Fc.statistics  / istats
    @views @. data.cΔt.data[:,itime] = data.cΔt.statistics / istats

    @. data.ϕ.statistics   = 0
    @. data.ϕf.statistics  = 0
    @. data.u.statistics   = 0
    @. data.v.statistics   = 0
    @. data.c.statistics   = 0
    @. data.Fc.statistics  = 0
    @. data.cΔt.statistics = 0

end

function savestats!(
    data :: ShallowWaterData2D,
    istats :: Int,
    itime :: Int,
)

    @views @. data.ϕ.data[:,:,itime]   = data.ϕ.statistics   / istats
    @views @. data.ϕf.data[:,:,itime]  = data.ϕf.statistics  / istats
    @views @. data.u.data[:,:,itime]   = data.u.statistics   / istats
    @views @. data.v.data[:,:,itime]   = data.v.statistics   / istats
    @views @. data.c.data[:,:,itime]   = data.c.statistics   / istats
    @views @. data.Fc.data[:,:,itime]  = data.Fc.statistics   / istats
    @views @. data.cΔt.data[:,:,itime] = data.cΔt.statistics / istats

    @. data.ϕ.statistics   = 0
    @. data.ϕf.statistics  = 0
    @. data.u.statistics   = 0
    @. data.v.statistics   = 0
    @. data.c.statistics   = 0
    @. data.Fc.statistics  = 0
    @. data.cΔt.statistics = 0

end

function savedata(
    data :: ShallowWaterData2D,
    model :: AbstractModel,
    dt    :: Real,
    itime :: Int,
    ifile :: Int,
    nsave :: Int
)

    fnc = joinpath(model.Directory,"day$(@sprintf("%05d",ifile)).nc")
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end
    if isfile(fnc)
        rm(fnc);
    end

    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) with YangShallowWater.jl",
    ))

    ds.dim["x"] = model.Grid.nx
    ds.dim["y"] = model.Grid.ny
    ds.dim["time"] = itime

    date = Date(0,1,1) + Day(ifile)

    ncx = defVar(ds,"x",Float32,("x",),attrib = Dict(
        "units"     => "meters",
        "long_name" => "x",
    ))

    ncy = defVar(ds,"y",Float32,("y",),attrib = Dict(
        "units"     => "meters",
        "long_name" => "y",
    ))

    nct = defVar(ds,"time",Int32,("time",),attrib = Dict(
        "units"     => "seconds since $(date) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    ncx[:] = model.Grid.x
    ncy[:] = model.Grid.y
    nct[:] = (dt * nsave) .* collect(1:itime)

    nc_ϕ = defVar(ds,"ϕ",Float32,("x","y","time"),attrib = Dict(
        "units"     => "m**2 s**-2",
        "long_name" => "geopotential",
    ))

    nc_ϕf = defVar(ds,"ϕf",Float32,("x","y","time"),attrib = Dict(
        "units"     => "m**2 s**-3",
        "long_name" => "large_scale_geopotential_tendency",
    ))

    nc_u = defVar(ds,"u",Float32,("x","y","time"),attrib = Dict(
        "units"     => "m s**-1",
        "long_name" => "x-direction wind",
    ))

    nc_v = defVar(ds,"v",Float32,("x","y","time"),attrib = Dict(
        "units"     => "m s**-1",
        "long_name" => "y-direction wind",
    ))

    nc_c = defVar(ds,"c",Float32,("x","y","time"),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "convection_grid_points",
    ))

    nc_Fc = defVar(ds,"Fc",Float32,("x","y","time"),attrib = Dict(
        "units"     => "m**2 s**-3",
        "long_name" => "geopotential_convection_tendency",
    ))

    nc_ct = defVar(ds,"cΔt",Float32,("x","y","time"),attrib = Dict(
        "units"     => "s",
        "long_name" => "time_left_for_convection",
    ))

    nc_ϕ[:]  = data.ϕ.data
    nc_ϕf[:] = data.ϕf.data
    nc_u[:]  = data.u.data
    nc_v[:]  = data.v.data
    nc_c[:]  = data.c.data
    nc_Fc[:] = data.Fc.data
    nc_ct[:] = data.cΔt.data

    close(ds)

end