module BNPMix

using Distributions
using ProgressMeter

import Base.run, Base.convert, StatsBase.sample

const Int   = Int64
const Float = Float64

include("utils.jl")
include("mcmc/metropolis.jl")
include("mcmc/sliceStepOut.jl")

include("xfamily/normalGammaIndependent.jl")
include("xfamily/normalNonConjugateHierarchySampled.jl")

include("mixture/cluster.jl")

include("nrmi/NGGP.jl")

include("mixture/mixture.jl")
include("mixture/mixtureNeal8.jl")
include("mixture/mixtureReuse.jl")
include("mixture/mixtureSlice.jl")

include("mcmc/sampler.jl")

end # module
