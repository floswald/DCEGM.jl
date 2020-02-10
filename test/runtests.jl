using DCEGM
using DelimitedFiles
using Test


@testset "Running DCEGM tests" begin

    include("basics.jl")
    include("linetests.jl")
    include("Envtests.jl")
    include("F_test.jl")

end
