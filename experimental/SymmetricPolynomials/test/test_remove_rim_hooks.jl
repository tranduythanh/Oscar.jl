using Test
using Oscar

# Include the SymmetricPolynomials module
include(joinpath(@__DIR__, "..", "src", "SymmetricPolynomials.jl"))
using .SymmetricPolynomials

@testset "Remove Rim Hooks" begin
    @testset "Test case 1" begin
        p = Partition([7,3])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([2,1])
    end

    @testset "Test case 2" begin
        p = Partition([7,2,1])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([1,1,1])
    end

    @testset "Test case 3" begin
        p = Partition([6,4])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([3])
    end

    @testset "Test case 4" begin
        p = Partition([6,3,1])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test isempty(result[1])
    end

    @testset "Test case 5" begin
        p = Partition([6,2,2])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([1,1,1])
    end

    @testset "Test case 6" begin
        p = Partition([5,4,1])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([3])
    end

    @testset "Test case 7" begin
        p = Partition([5,3,2])
        result = remove_rim_hooks(p, 7, (3, 4))
        @test result[1] == Partition([2,1])
    end

    @testset "Test case 8" begin
        p = Partition([5,5])
        result = remove_rim_hooks(p, 5, (2, 3))
        @test isempty(result[1])
    end

    @testset "Additional properties" begin
        for test_case in [
            ([7,3], 7, (3, 4)),
            ([7,2,1], 7, (3, 4)),
            ([6,4], 7, (3, 4)),
            ([6,3,1], 7, (3, 4)),
            ([6,2,2], 7, (3, 4)),
            ([5,4,1], 7, (3, 4)),
            ([5,3,2], 7, (3, 4)),
            ([5,5], 5, (2, 3))
        ]
            p, rim_size, grid = test_case
            result = remove_rim_hooks(Partition(p), rim_size, grid)
            
            # Check that the resulting partition fits within the grid
            @test isempty(result[1]) || length(result[1]) <= grid[1]
            @test isempty(result[1]) || result[1][1] <= grid[2]
            
            # Check that the number of removed hooks is non-negative
            @test result[2] >= 0
            
            # Check that the skew partition height is non-negative
            @test result[3] >= 0
        end
    end
end
