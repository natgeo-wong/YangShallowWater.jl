function calcYangLargeScale(ϕ, ϕforcing, clock)

    ϕf = ϕforcing.ϕl
    if ϕforcing.wtg
        ϕf += (mean(ϕ) - ϕforcing.ϕ0) / ϕforcing.τl
    end

    return ϕf * clock.dt

end