using YangShallowWater
using Test

@testset "YangShallowWater.jl" begin
    nx, Lx = 10, 2.0
    grid = GenerateGrid(; nx, Lx)
    @test_throws "type OneDGrid has no field a" grid.Krsq

    nx, Lx = 8, 2.0
    ny, Ly = 12, 3.0
    grid = GenerateGrid(; nx, ny, Lx, Ly)
    @test size(grid.Ksq) == (nx, ny)
    @test size(grid.Krsq) == (Int(nx/2+1), ny)
end
