

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
    @devzeros Dev Complex{T} (grid.nk, grid.nl) uh vh ϕh Fch

    return Vars(u, v, ϕ, Fc, uh, vh, ηh, Fch)

end

function calcN!(N, sol, t, clock, vars, params, grid)

    @. vars.uh  = sol[:,:,1]
    @. vars.vh  = sol[:,:,2]
    @. vars.ϕh  = sol[:,:,3]
    @. vars.Fch = sol[:,:,4]

    @. N[:, 1] = - im * grid.nk * vars.ϕh - vars.uh / τd  # - ∂ϕ/∂x - u/τd
    @. N[:, 2] = - im * grid.nl * vars.ϕh - vars.uh / τd  # - ∂ϕ/∂y - v/τd
    @. N[:, 3] = - im * (grid.nk * vars.uh + grid.nl * vars.vh) + vars.Fch

    dealias!(N, grid)

    return nothing

end

function Equation(dev, params, grid)

    T = eltype(grid)
    L = zeros(dev, T, (grid.nkr, 3))
    D = @. - params.ν * grid.kr^(2*params.nν)

    L[:, 1] .= D # for u equation
    L[:, 2] .= D # for v equation
    L[:, 3] .= D # for η equation

    return FourierFlows.Equation(L, calcN!, grid)

end

function updatevars!(prob)
    
    vars, grid, sol = prob.vars, prob.grid, prob.sol

    @. vars.uh = sol[:,:,1]
    @. vars.vh = sol[:,:,2]
    @. vars.ηh = sol[:,:,3]

    ldiv!(vars.u, grid.rfftplan, deepcopy(sol[:, 1])) # use deepcopy() because irfft destroys its input
    ldiv!(vars.v, grid.rfftplan, deepcopy(sol[:, 2])) # use deepcopy() because irfft destroys its input
    ldiv!(vars.η, grid.rfftplan, deepcopy(sol[:, 3])) # use deepcopy() because irfft destroys its input

    return nothing

end

function set_uvη!(prob, u0, v0, η0)

    vars, grid, sol = prob.vars, prob.grid, prob.sol

    A = typeof(vars.u) # determine the type of vars.u

    mul!(vars.uh, grid.rfftplan, A(u0))
    mul!(vars.vh, grid.rfftplan, A(v0))
    mul!(vars.ηh, grid.rfftplan, A(η0))

    @. sol[:, 1] = vars.uh
    @. sol[:, 2] = vars.vh
    @. sol[:, 3] = vars.ηh

    updatevars!(prob)

  return nothing

end