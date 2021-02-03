using CooperationPolarization
using Test

@testset "CooperationPolarization.jl" begin

    # setup.jl: random_sets
    @test random_sets(2, 1, 1, 2).p == [0, 1]
    @test random_sets(8, 3, 2, 2).p == [0, 0, 0, 0, 1, 1, 1, 1]
    @test random_sets(8, 3, 2).p == random_sets(8, 3, 2, 2).p 

    # setup.jl: vectobase3
    @test vectobase3([0; 0; 0]) == 0
    @test vectobase3([0; 1; 1]) == 4
    @test vectobase3([0;-1;-1]) == -4

end
