using NRMMM
using Base.Test

@testset "metropolis"          begin include("metropolisTest.jl")          end
@testset "sliceFiniteInterval" begin include("sliceFiniteIntervalTest.jl") end
@testset "sliceStepOut"        begin include("sliceStepOutTest.jl")        end
