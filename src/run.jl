function run(
    model :: Model2DSpectral;
    dt :: Real = 5.,
    nsteps :: Int,
    nstats :: Int = Int(60 /dt),
    nsave  :: Int = Int(600/dt),
    ϕ0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
    u0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
    v0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
)


    if !iszero(mod(30,dt))
        error("$(modulelog()) - timestep dt=$nstats should be a factor of 30 (30 seconds is the maximum timestep currently allowed")
    end

    E = Equation(model)
    prob = FourierFlows.Problem(
        E, "FilteredRK4", dt, 
        model.Grid, model.Variables, model.Parameters
    )
    set_uvϕ!(prob, u0, v0, ϕ0, model.Grid, model.Variables)

    ntotal, data = createsavedata(dt,nsteps,nstats,nsave,model.Grid)

    is = 1
    istats = 0
    itime  = 0
    ifile  = 0
    @showprogress "Running YangShallowWater.jl over $nsteps steps, dt = $dt s, model elapsed time = $(nsteps*dt/86400) days ..." for it = 1 : nsteps
        stepforward!(prob)

        if iszero(mod(it,nstats))
            istats += 1
            updatevars!(prob,model.Grid,model.Variables)
            updatestats!(data,prob,prob.params)
            if iszero(mod(it,nsave))
                itime += 1
                savestats!(data,istats,itime)
                istats = 0
                if iszero(mod(itime,ntotal))
                    ifile += 1
                    savedata(data,model,dt,itime,ifile,nsave)
                    itime = 0
                end
            end
        end

    end

    if !iszero(itime)
        ifile += 1
        savedata(data,model,dt,itime,ifile,nsave)
    end

    return 

end