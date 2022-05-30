using BrillNoetherGeneralLimitSpinBundleCounter
using Test

A = matrix(ZZ,[[64,0],[0,0],[0,512],[0,1024],[0,896],[0,1024],[0,512],[0,0],[0,64]])

@testset "Compute some Brill-Noether numbers" begin
    @test Counter([[1,1]]) == matrix(ZZ, [[2,0],[0,2]])
    @test Counter([[1,3]]) == matrix(ZZ, [[1],[0]])
    @test Counter([[1,2],[1,2]]) == matrix(ZZ,[[2,0],[0,0],[0,2]])
    @test Counter([[1,2],[1,2],[1,2],[1,2],[2,3],[2,3],[1,3],[1,3]]) == A
    @test Generic_Brill_Noether_numbers([[1,2],[1,2],[1,2],[1,2],[2,3],[2,3],[1,3],[1,3]]) == [64,4032]
end
