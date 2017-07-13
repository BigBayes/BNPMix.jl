module NRMMM

using Distributions
using ProgressMeter

import Base.run, Base.convert, StatsBase.sample

const Int   = Int64
const Float = Float64

#include("nrmi.jl")
include("mcmc/metropolis.jl")
include("mcmc/sliceStepOut.jl")
include("mcmc/ess.jl")

include("xfamily/normalGammaIndependent.jl")
include("xfamily/normalNonConjugateHierarchySampled.jl")

include("mixture/cluster.jl")

include("nrmi/NGGP.jl")

include("mixture/mixture.jl")
include("mixture/mixtureNeal8.jl")
include("mixture/mixtureReuse.jl")

include("mcmc/sampler.jl")

end # module
