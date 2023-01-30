function calcYangConvection!(ch, c, ϕ, convection, clock, grid, ::SpectralVars)

    calcYangConvection!(c, ϕ, convection, clock, grid)

    mul!(ch, grid.rfftplan, c)

    return nothing

end

function calcYangConvection!(c, ϕ, convection, clock, grid)

    updateYangConvectionLocation!(convection, ϕ)
    updateYangConvectionFlux!(c, convection, clock, grid)

    return nothing

end

function updateYangConvectionLocation!(convection, ϕ)

    old = convection.ConvectionFlux.old
    new = convection.ConvectionFlux.new
    Δt  = convection.ConvectionFlux.Δt
    
    for ii in eachindex(old)
        if isone(old[ii])
            Δt[ii] -= clock.dt
            if Δt[ii] < 0
                Δt[ii] = 0
                old[ii] = 0
            end
        end
    end

    @. new = ϕ > convection.ϕc
    
    for ii in eachindex(new)
        if isone(new[ii] - old[ii])
            Δt[ii] = params.convection.τc
        end
    end
    
    return nothing

end

function updateYangConvectionFlux!(c, convection::Convection2D, clock, G::TwoDGrid)

    c_new  = convection.ConvectionFlux.new
    c_Δt   = convection.ConvectionFlux.Δt

    CG.nx = convection.ConvectionGrid.nx
    CG.ny = convection.ConvectionGrid.ny
    CG.cx = convection.ConvectionGrid.cx
    CG.cy = convection.ConvectionGrid.cy
    CG.Fc = convection.ConvectionGrid.Fc

    τcd2 = convection.τc/2

    for iy = 1 : G.ny, ix = 1 : G.nx

        if isone(c_new[ix,iy])

            τ = 1 - (c_Δt[ix,iy]/τcd2-1)^2

            for iiy = 1 : CG.ny, iix = 1 : CG.nx
                ic_x = Int(mod(ix+CG.cx[iix],G.nx)); if iszero(ic_x); ic_x = G.nx end
                ic_y = Int(mod(iy+CG.cy[iiy],G.ny)); if iszero(ic_y); ic_y = G.ny end
                c[ic_x,ic_y] = CG.Fc[iix,iiy] * τ * clock.dt
            end

        end

    end

    return nothing

end

function updateYangConvectionFlux!(c, convection::Convection1D, clock, G::OneDGrid)

    c_new  = convection.ConvectionFlux.new
    c_Δt   = convection.ConvectionFlux.Δt

    CG.nx = convection.ConvectionGrid.nx
    CG.cx = convection.ConvectionGrid.cx
    CG.Fc = convection.ConvectionGrid.Fc

    τcd2 = convection.τc/2

    for ix = 1 : G.nx

        if isone(c_new[ix])

            τ = 1 - (c_Δt[ix]/τcd2-1)^2

            for iix = 1 : CG.nx
                ic_x = Int(mod(ix+CG.cx[iix],G.nx)); if iszero(ic_x); ic_x = G.nx end
                c[ic_x] = CG.Fc[iix] * τ * clock.dt
            end

        end

    end

    return nothing

end