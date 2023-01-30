function run(
    model :: Model2DSpectral;
    dt :: Real = 5.,
    nsteps :: Real,
    ϕ0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
    u0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
    v0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
) where FT <: Real

    E = Equation(model)
    prob = FourierFlows.Problem(
        E, "FilteredRK4", dt, 
        model.Grid, model.Variables, model.Parameters
    )
    set_uvϕ!(prob, u0, v0, ϕ0, model.Grid, model.Variables)

    for it = 0 : nsteps
        stepforward!(prob)
    end
    updatevars!(prob, model.Grid, model.Variables)

    return prob

end