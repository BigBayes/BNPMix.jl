using Distributions
include("normalGammaIndependent.jl")

mutable struct NormalNonConjugateHierarchySampled <: Hierarchy
  prior::NormalGammaIndependent
  sumX::Float64
  sumXX::Float64
  number::Int64
  param::Distributions.Normal{Float64}
end

NormalNonConjugateHierarchySampled(prior::NormalGammaIndependent) =
NormalNonConjugateHierarchySampled(prior, 0.0, 0.0, 0, drawSample(prior))

function addDatum(h::NormalNonConjugateHierarchySampled, datum::Float64)
  h.number += 1
	h.sumX += datum
	h.sumXX += datum*datum
end

function removeDatum(h::NormalNonConjugateHierarchySampled, datum::Float64)
  h.number -= 1
	h.sumX -= datum
	h.sumXX -= datum*datum
  assert(h.number >= 0)
  assert(h.sumXX >= -1e-10)
end

function sample(h::NormalNonConjugateHierarchySampled)
  #print("NormalNonConjugateHierarchySampled - sample")
  if (h.number == 0)
    h.param = drawSample(h.prior)
  else
    paramMean = mean(h.param)
    paramPrecision = 1.0/var(h.param)
    newPrecision = h.prior.meanPrecision + h.number*paramPrecision
    newMean = (h.prior.meanPrecision*h.prior.meanMean + paramPrecision*h.sumX) / newPrecision
    newParamMean = rand(Normal(newMean, 1.0/sqrt(newPrecision)))
    newShape = h.prior.precisionShape +.5*h.number
    newInvScale = h.prior.precisionInvScale +
              .5*(h.sumXX - 2.0*paramMean*h.sumX + h.number*paramMean*paramMean)
    newParamPrecision = rand(Gamma(newShape, 1.0/newInvScale))
    h.param = Normal(newParamMean, 1.0/sqrt(newParamPrecision))
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

function logPredictive(h::NormalNonConjugateHierarchySampled, datum::Float64)
  return logpdf(h.param, datum)
end

######### NormalNonConjugateFactorySampled

abstract type Factory end

struct NormalNonConjugateFactorySampled <: Factory
end

function construct{P<:Prior}(f::NormalNonConjugateFactorySampled, prior::P)
  return NormalNonConjugateHierarchySampled(prior)
end
