function run(
    P :: YSWParams{FT},
    G :: TwoDGrid;
    dt :: Real = 5.,
    nsteps :: Real,
    ϕ0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
    u0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
    v0  :: Array{<:Real,2} = zeros(FT,G.nx,G.ny),
) where FT <: Real

    E = Equation(P,G)
    V = Vars(G)
    prob = FourierFlows.Problem(E, "FilteredRK4", dt, G, V, P)
    set_uvϕ!(prob, u0, v0, ϕ0)

    for it = 0 : nsteps
        stepforward!(prob)
    end
    updatevars!(prob)

    return prob

end