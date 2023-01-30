function calcN!(N, sol, t, clock, vars :: SpectralVars, params :: SimpleParams , grid :: TwoDGrid)
    
    dealias!(sol, grid)

    vars.uh  .= view(sol, :, :, 1)
    vars.vh  .= view(sol, :, :, 2)
    vars.ϕh  .= view(sol, :, :, 3)

    @views @. N[:, :, 1] = - im * grid.kr * vars.ϕh                                    # - ∂ϕ/∂x
    @views @. N[:, :, 2] = - im * grid.l  * vars.ϕh                                    # - ∂ϕ/∂y
    @views @. N[:, :, 3] = - im * (grid.kr * vars.uh + grid.l * vars.vh) * params.c^2  # - c^2 * (∂u/∂x + ∂v/∂y)

    return nothing

end

function calcN!(N, sol, t, clock, vars :: SpectralVars, params :: ForcingParams, grid :: TwoDGrid)
    
    dealias!(sol, grid)

    vars.uh  .= view(sol, :, :, 1)
    vars.vh  .= view(sol, :, :, 2)
    vars.ϕh  .= view(sol, :, :, 3)

    @views @. N[:, :, 1] = - im * grid.kr * vars.ϕh                                    # - ∂ϕ/∂x
    @views @. N[:, :, 2] = - im * grid.l  * vars.ϕh                                    # - ∂ϕ/∂y
    @views @. N[:, :, 3] = - im * (grid.kr * vars.uh + grid.l * vars.vh) * params.c^2  # - c^2 * (∂u/∂x + ∂v/∂y)

    addforcing!(N, sol, t, clock, vars, params, grid)

    return nothing

end

function addforcing!(N, sol, t, clock, vars :: SpectralVars, params, grid :: TwoDGrid)

    ldiv!(vars.ϕ, grid.rfftplan, deepcopy(sol[:, :, 3]))

    if params.forcing
        ϕf = params.ϕforcing.Equation(vars.ϕ, params.ϕforcing, clock)
    end
    
    # call calcF! to compute ch and store it in vars.ch
    if params.convection
        params.convection.Equation(
            vars.ch, vars.c, vars.ϕ,
            params.convection, clock, grid, vars
        )
    end
    
    # add ch on the nonlinear term for ϕ
    @views @. N[:, :, 3] += vars.ch + ϕf
  
    return nothing
end

function Equation(model::Model2DSpectral)

    grid = model.Grid
    params = model.Parameters

    T = eltype(grid)
    dev = grid.device

    L = zeros(dev, T, (grid.nkr, grid.nl, 3))

    D = @. - params.ν * grid.Krsq^params.nν  # - ν (k²+l²)ⁿ

    @views @. L[:, :, 1] = D - 1 / params.τd    # for u
    @views @. L[:, :, 2] = D - 1 / params.τd    # for v

    return FourierFlows.Equation(L, calcN!, grid)

end

function updatevars!(prob, ::TwoDGrid, ::SpectralVars)
    
    vars, grid, sol = prob.vars, prob.grid, prob.sol

    vars.uh  .= sol[:, :, 1]
    vars.vh  .= sol[:, :, 2]
    vars.ϕh  .= sol[:, :, 3]

    # use deepcopy() below because irfft destroys its input
    ldiv!(vars.u, grid.rfftplan, deepcopy(sol[:, :, 1]))
    ldiv!(vars.v, grid.rfftplan, deepcopy(sol[:, :, 2]))
    ldiv!(vars.ϕ, grid.rfftplan, deepcopy(sol[:, :, 3]))
    ldiv!(vars.c, grid.rfftplan, deepcopy(vars.ch))

    return nothing

end

function set_uvϕ!(prob, u0, v0, ϕ0, G::TwoDGrid, V::SpectralVars)

    vars, grid, sol = prob.vars, prob.grid, prob.sol

    A = typeof(vars.u) # determine the type of vars.u

    mul!(vars.uh, grid.rfftplan, A(u0))
    mul!(vars.vh, grid.rfftplan, A(v0))
    mul!(vars.ϕh, grid.rfftplan, A(ϕ0))

    @views sol[:, :, 1] .= vars.uh
    @views sol[:, :, 2] .= vars.vh
    @views sol[:, :, 3] .= vars.ϕh

    updatevars!(prob, G, V)

  return nothing

end
