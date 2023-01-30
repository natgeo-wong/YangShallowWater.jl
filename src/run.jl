function run(
    model :: Model2DSpectral;
    dt :: Real = 5.,
    nsteps :: Int,
    nsave  :: Int,
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
    end

    updatevars!(prob, model.Grid, model.Variables)
    ϕf[:,:,end] .= prob.vars.ϕ
    cf[:,:,end] .= prob.vars.c
    tf[:,:,end] .= prob.params.convection.ConvectionFlux.Δt

    return ϕf,cf,tf,prob.clock

end