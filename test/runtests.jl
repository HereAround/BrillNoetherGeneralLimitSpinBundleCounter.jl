using BrillNoetherGeneralLimitSpinBundleCounter
using Test

@testset "Compute some Brill-Noether numbers" begin
    @test Counter([[1,2],[1,2]]) == matrix(ZZ,[[2,0],[0,0],[0,2],[2,2]])
end
