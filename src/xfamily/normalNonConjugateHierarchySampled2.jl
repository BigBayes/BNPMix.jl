using Distributions

export
    NormalNonConjugateFactorySampled,
    construct

mutable struct NormalNonConjugateHierarchySampled <: Hierarchy
  prior  ::  NormalGammaIndependent # Gamma-Independent prior on the Normal distribution
  sumX   ::  Float                  # Running mean of data belonging to the associated cluster
  sumXX  ::  Float                  # Running variance of data belonging to the associated cluster
  number ::  Int                    # Number of data belonging to the associated cluster
  param  ::  Distributions.Normal{Float}
end

NormalNonConjugateHierarchySampled(prior::NormalGammaIndependent) =
NormalNonConjugateHierarchySampled(prior, 0.0, 0.0, 0, drawSample(prior))

function addDatum(h::NormalNonConjugateHierarchySampled, datum::Array{Float, 1})
  h.number += 1
	h.sumX += datum[1]
	h.sumXX += datum[1]*datum[1]
end

function removeDatum(h::NormalNonConjugateHierarchySampled, datum::Array{Float, 1})
  h.number -= 1
	h.sumX -= datum[1]
	h.sumXX -= datum[1]*datum[1]
  assert(h.number >= 0)
  assert(h.sumXX >= -1e-10)
end

function sample(h::NormalNonConjugateHierarchySampled)
  if (h.number == 0)
    h.param = drawSample(h.prior)
  else
    paramMean = mean(h.param)
    paramPrecision = 1.0/var(h.param)
    newPrecision = h.prior.meanPrecision + h.number*paramPrecision
    newMean = (h.prior.meanPrecision*h.prior.meanMean + paramPrecision*h.sumX) / newPrecision
    newParamMean = rand(Normal(newMean, 1.0/sqrt(newPrecision)))
    h.param = Normal(newParamMean, 1.0/sqrt(paramPrecision))
  end
end

function logJoint(h::NormalNonConjugateHierarchySampled)
  neghalflog2pi = -.5*log(2*pi)
  paramprecision = 1.0/var(h.param)
  parammean = mean(h.param)
  return h.number * (neghalflog2pi + .5*log(paramprecision)) -
    .5*paramprecision*(h.sumXX - 2.0*h.sumX*parammean + h.number*parammean*parammean) +
    logProbability(h.prior, h.param)
end

function logPredictive(h::NormalNonConjugateHierarchySampled, datum::Array{Float, 1})
  return logpdf(h.param, datum[1])
end

######### NormalNonConjugateFactorySampled

abstract type Factory end

struct NormalNonConjugateFactorySampled <: Factory
end

function construct{P<:Prior}(f::NormalNonConjugateFactorySampled, prior::P)
  return NormalNonConjugateHierarchySampled(prior)
end
