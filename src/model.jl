
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

function Vars(::Dev, grid) where Dev

    T = eltype(grid)
    @devzeros Dev T (grid.nx, grid.ny) u v ϕ Fc
    @devzeros Dev Complex{T} (grid.nkr, grid.nl) uh vh ϕh Fch

    return Vars(u, v, ϕ, Fc, uh, vh, ηh, Fch)

end

function calcN!(N, sol, t, clock, vars, params, grid)

    @. vars.uh  = view(sol, :, :, 1)
    @. vars.vh  = view(sol, :, :, 2)
    @. vars.ϕh  = view(sol, :, :, 3)
    @. vars.Fch = view(sol, :, :, 4)

    @views @. N[:, :, 1] = - im * grid.kr * vars.ϕh                                    # - ∂ϕ/∂x
    @views @. N[:, :, 2] = - im * grid.l  * vars.ϕh                                    # - ∂ϕ/∂y
    @views @. N[:, :, 3] = - im * (grid.kr * vars.uh + grid.l * vars.vh) + vars.Fch    # - ∂u/∂x - ∂v/∂y - Fc
    @views @. N[:, :, 4] = 0                                                           # 0 for now

    dealias!(N, grid)

    return nothing

end

function Equation(dev, params, grid)

    T = eltype(grid)
    L = zeros(dev, T, (grid.nkr, grid.l, 3))
    D = @. - params.ν * grid.Krsq^params.nν

    # for u, v, φ
    for j in 1:3
        L[:, :, j] .= D
    end

    @. L[:, :, 1] = - 1 / τd    # the -1/τd drag on u
    @. L[:, :, 2] = - 1 / τd    # the -1/τd drag on v

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