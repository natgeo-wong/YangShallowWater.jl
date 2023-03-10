function GenerateGrid(
	FT = Float64;
	nx :: Int,
    ny :: Int = 1,
    Lx :: Real,
    Ly :: Real = 0.
)

    if isone(ny) || iszero(Ly)
        return OneDGrid(;
            nx, Lx, x0 = Lx/(2nx),
            T = FT
        )
    else
        return TwoDGrid(;
            nx, Lx, x0 = Lx/(2nx),
            ny, Ly, y0 = Ly/(2ny),
            T = FT
        )
    end

end