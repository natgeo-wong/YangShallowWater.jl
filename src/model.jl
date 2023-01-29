struct Vars{Aphys, Atrans} <: AbstractVars
      u :: Aphys
      v :: Aphys
      ϕ :: Aphys
     Fc :: Aphys
     uh :: Atrans
     vh :: Atrans
     ϕh :: Atrans
    Fch :: Atrans
end

function Vars(grid)

    T = eltype(grid)
    Dev = typeof(grid.device)

    @devzeros Dev T (grid.nx, grid.ny) u v ϕ Fc
    @devzeros Dev Complex{T} (grid.nkr, grid.nl) uh vh ϕh Fch

    return Vars(u, v, ϕ, Fc, uh, vh, ϕh, Fch)

end

function calcN!(N, sol, t, clock, vars, params, grid)

    @. vars.uh  = view(sol, :, :, 1)
    @. vars.vh  = view(sol, :, :, 2)
    @. vars.ϕh  = view(sol, :, :, 3)
    @. vars.Fch = view(sol, :, :, 4)

    @views @. N[:, :, 1] = - im * grid.kr * vars.ϕh                                    # - ∂ϕ/∂x
    @views @. N[:, :, 2] = - im * grid.l  * vars.ϕh                                    # - ∂ϕ/∂y
    @views @. N[:, :, 3] = - im * (grid.kr * vars.uh + grid.l * vars.vh) * params.c^2 - vars.Fch    # - c^2 * (∂u/∂x - ∂v/∂y) - Fc
    @views @. N[:, :, 4] = 0                                                           # 0 for now

    dealias!(N, grid)

    return nothing

end

function Equation(params, grid)

    T = eltype(grid)
    dev = grid.device

    L = zeros(dev, T, (grid.nkr, grid.nl, 3))

    D = @. - params.ν * grid.Krsq^params.nν  # - ν (k²+l²)ⁿ

    @views @. L[:, :, 1] = D - 1 / params.τd    # for u
    @views @. L[:, :, 2] = D - 1 / params.τd    # for v
    @views @. L[:, :, 3] = D                    # for ϕ

    return FourierFlows.Equation(L, calcN!, grid)

end

function updatevars!(prob)
    
    vars, grid, sol = prob.vars, prob.grid, prob.sol

    @. vars.uh = sol[:, :, 1]
    @. vars.vh = sol[:, :, 2]
    @. vars.ϕh = sol[:, :, 3]
    @. vars.Fch = sol[:, :, 4]

    # use deepcopy() below because irfft destroys its input
    ldiv!(vars.u, grid.rfftplan, deepcopy(sol[:, :, 1]))
    ldiv!(vars.v, grid.rfftplan, deepcopy(sol[:, :, 2]))
    ldiv!(vars.ϕ, grid.rfftplan, deepcopy(sol[:, :, 3]))
    ldiv!(vars.Fc, grid.rfftplan, deepcopy(sol[:, :, 4]))

    return nothing

end

function set_uvϕFc!(prob, u0, v0, ϕ0, Fc0)

    vars, grid, sol = prob.vars, prob.grid, prob.sol

    A = typeof(vars.u) # determine the type of vars.u

    mul!(vars.uh, grid.rfftplan, A(u0))
    mul!(vars.vh, grid.rfftplan, A(v0))
    mul!(vars.ϕh, grid.rfftplan, A(ϕ0))
    mul!(vars.Fch, grid.rfftplan, A(Fc0))

    @views @. sol[:, :, 1] = vars.uh
    @views @. sol[:, :, 2] = vars.vh
    @views @. sol[:, :, 3] = vars.ϕh
    @views @. sol[:, :, 4] = vars.Fh

    updatevars!(prob)

  return nothing

end