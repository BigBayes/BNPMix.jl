using NRMMM
using Base.Test

@testset "metropolis"          begin include("metropolisTest.jl")          end
@testset "sliceFiniteInt64erval" begin include("sliceFiniteInt64ervalTest.jl") end
@testset "sliceStepOut"        begin include("sliceStepOutTest.jl")        end
