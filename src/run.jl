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

    ϕf = zeros(model.Grid.nx,model.Grid.ny,Int(ceil(nsteps/nsave)+1))
    cf = zeros(model.Grid.nx,model.Grid.ny,Int(ceil(nsteps/nsave)+1))
    tf = zeros(model.Grid.nx,model.Grid.ny,Int(ceil(nsteps/nsave)+1))
    ϕf[:,:,1] .= ϕ0

    if !iszero(mod(nsave,nstats))
        error("$(modulelog()) - nstats=$nstats is not a factor of nsave=$nsave")
    end

    is = 1
    for it = 1 : nsteps
        stepforward!(prob)
        if iszero(mod(it,nsave))
            updatevars!(prob, model.Grid, model.Variables)
            is += 1
            ϕf[:,:,is] .= prob.vars.ϕ
            cf[:,:,is] .= prob.vars.c
            tf[:,:,is] .= prob.params.convection.ConvectionFlux.Δt
        end
        if iszero(mod(it,nlogs))
            @info "$(modulelog()) - Step $it of $nsteps"
        end
    end

    @info "$(modulelog()) - Step $nsteps of $nsteps"
    updatevars!(prob, model.Grid, model.Variables)
    ϕf[:,:,end] .= prob.vars.ϕ
    cf[:,:,end] .= prob.vars.c
    tf[:,:,end] .= prob.params.convection.ConvectionFlux.Δt

    return ϕf,cf,tf,prob.clock

end