using CooperationPolarization
using Test

@testset "CooperationPolarization.jl" begin

    # setup.jl: random_sets
    @test random_sets(2, 1, 1, 2).affiliations == [1, 2]
    @test random_sets(8, 3, 2, 2).affiliations == [1, 1, 1, 1, 2, 2, 2, 2]
    @test random_sets(8, 3, 2).affiliations == random_sets(8, 3, 2, 2).affiliations

    # setup.jl: vectobase3
    @test vectobase3([0; 0; 0]) == 0
    @test vectobase3([0; 1; 1]) == 4
    @test vectobase3([0;-1;-1]) == -4

end
