function run(
    model :: Model2DSpectral;
    dt :: Real = 5.,
    nsteps :: Int,
    nstats :: Int,
    nsave  :: Int,
    nlogs  :: Int = 100,
    ϕ0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
    u0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
    v0  :: Array{<:Real,2} = zeros(model.Grid.nx,model.Grid.ny),
)

    E = Equation(model)
    prob = FourierFlows.Problem(
        E, "FilteredRK4", dt, 
        model.Grid, model.Variables, model.Parameters
    )
    set_uvϕ!(prob, u0, v0, ϕ0, model.Grid, model.Variables)

    if !iszero(mod(nsave,nstats))
        error("$(modulelog()) - nstats=$nstats is not a factor of nsave=$nsave")
    end

    is = 1
    @showprogress "Running YangShallowWater.jl over $nsteps steps, dt = $dt s, model elapsed time = $(nsteps*dt/86400) days ..." for it = 1 : nsteps
        stepforward!(prob)
    end

    return 

end