function calcYangConvection!(ch, c, ϕ, ϕf, convection, clock, grid, ::SpectralVars)

    calcYangConvection!(c, ϕ, ϕf, convection, clock, grid)

    mul!(ch, grid.rfftplan, c)

    return nothing

end

function calcYangConvection!(c, ϕ, ϕf, convection, clock, grid)

    updateYangConvectionLocation!(convection, clock, ϕ)
    updateYangConvectionFlux!(c, ϕf, convection, clock, grid)

    return nothing

end

function updateYangConvectionLocation!(convection, clock, ϕ)

    old = convection.ConvectionFlux.old
    new = convection.ConvectionFlux.new
    Δt  = convection.ConvectionFlux.Δt

    @. old = new
    
    for ii in eachindex(Δt)
        if Δt[ii] < 0
            Δt[ii] = 0
        else
            Δt[ii] -= clock.dt
        end
    end

    @. new = ϕ > convection.ϕc
    
    for ii in eachindex(new)
        if (isone(new[ii] - old[ii])) && iszero(Δt[ii]) 
            Δt[ii] = convection.τc
        end
    end
    
    return nothing

end

function updateYangConvectionFlux!(c, ϕf, convection::Convection2D, clock, G::TwoDGrid)

    c .= 0
    c_new  = convection.ConvectionFlux.new
    c_Δt   = convection.ConvectionFlux.Δt

    CG = convection.ConvectionGrid

    τcd2 = convection.τc/2
    q = ϕf./convection.Sc

    for iy = 1 : G.ny, ix = 1 : G.nx

        if !iszero(c_Δt[ix,iy])

            τ = (1 - (c_Δt[ix,iy]/τcd2-1)^2) * q

            for iiy = 1 : CG.ny, iix = 1 : CG.nx
                ic_x = Int(mod(ix+CG.cx[iix],G.nx)); if iszero(ic_x); ic_x = G.nx end
                ic_y = Int(mod(iy+CG.cy[iiy],G.ny)); if iszero(ic_y); ic_y = G.ny end
                c[ic_x,ic_y] += CG.Fc[iix,iiy] * τ
            end

        end

    end

    return nothing

end

function updateYangConvectionFlux!(c, ϕf, convection::Convection1D, clock, G::OneDGrid)

    c .= 0
    c_new  = convection.ConvectionFlux.new
    c_Δt   = convection.ConvectionFlux.Δt

    CG = convection.ConvectionGrid

    τcd2 = convection.τc/2
    q = ϕf./convection.Sc

    for ix = 1 : G.nx

        if isone(c_new[ix])

            τ = (1 - (c_Δt[ix]/τcd2-1)^2) * q

            for iix = 1 : CG.nx
                ic_x = Int(mod(ix+CG.cx[iix],G.nx)); if iszero(ic_x); ic_x = G.nx end
                c[ic_x] += CG.Fc[iix] * τ
            end

        end

    end

    return nothing

end