using DCEGM
using Test
# using TestSetExtensions


@testset "Running DCEGM tests" begin

    # @includetests ARGS
    include("linetests.jl")
    include("Envtests.jl")

end
    


