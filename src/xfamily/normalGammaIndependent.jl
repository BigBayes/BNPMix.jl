using Distributions
include("../mcmc/metropolis.jl")

abstract type Hierarchy end
abstract type Prior end

mutable struct NormalGammaIndependent <: Prior
    meanMean::Float64
    meanPrecision::Float64
    precisionShape::Float64
    precisionInvScaleAlpha::Float64
    precisionInvScaleBeta::Float64
    precisionInvScale::Float64
    data::Union{Set{Hierarchy}, Void}

    function NormalGammaIndependent(meanMean::Float64, meanPrecision::Float64, precisionShape::Float64, precisionInvScaleAlpha::Float64, precisionInvScaleBeta::Float64, precisionInvScale::Float64)
    	assert((meanMean >= zero(meanMean)) && (meanPrecision >= zero(meanPrecision)) && (precisionShape >= zero(precisionShape)) &&
      (precisionInvScaleAlpha >= zero(precisionInvScaleAlpha)) && (precisionInvScaleBeta >= zero(precisionInvScaleBeta)) &&
      (precisionInvScale >= zero(precisionInvScale)))
    	new(meanMean, meanPrecision, precisionShape, precisionInvScaleAlpha, precisionInvScaleBeta, precisionInvScale, nothing)
    end
end

NormalGammaIndependent(meanMean::Float64, meanPrecision::Float64, precisionShape::Float64, precisionInvScaleAlpha::Float64, precisionInvScaleBeta::Float64) =
NormalGammaIndependent(meanMean, meanPrecision, precisionShape, precisionInvScaleAlpha, precisionInvScaleBeta, precisionInvScaleAlpha/precisionInvScaleBeta)

function drawSample(d::NormalGammaIndependent)
  return Normal(rand(Normal(d.meanMean, 1.0/sqrt(d.meanPrecision))), 1.0/sqrt(rand(Gamma(d.precisionShape, 1.0/d.precisionInvScale))) )
end

function logNormalizer(d::NormalGammaIndependent)
  halfLogTwoPi = .5*log(2*pi)
  return halfLogTwoPi -.5*log(d.meanPrecision) -
            d.precisionShape*log(d.precisionInvScale) +
            lgamma(d.precisionShape)
end

function logProbability(d::NormalGammaIndependent, normal::Distributions.Normal{Float64})
  diff = mean(normal) - d.meanMean
  normalprecision = 1.0/var(normal)
	return  -.5*d.meanPrecision*diff*diff +
          (d.precisionShape-1)*log(normalprecision) -
          d.precisionInvScale*normalprecision -
					logNormalizer(d)
end

### Sampling

function sample(d::NormalGammaIndependent, data::Set{Hierarchy})
  #print("NormalGammaIndependent - sample\n")
  d.data = data
  #sampleMeanMean(d)
  #sampleMeanPrecision(d)
  #samplePrecisionShape(d)
  samplePrecisionInvScale(d)
end

function samplePrecisionInvScale(d::NormalGammaIndependent)

  function propose(precisionInvScale::Float64)
    return precisionInvScale * exp(.5 * rand(Normal()))
  end

  function logratio(isnew::Float64, isold::Float64)
    d.precisionInvScale = isold
    result = d.precisionInvScaleAlpha * (log(isnew) - log(isold)) -
    d.precisionInvScaleBeta * (isnew - isold)
    for datum in d.data
      result -= logJoint(datum)
    end
    d.precisionInvScale = isnew
    for datum in d.data
      result += logJoint(datum)
    end
    return result
  end

  sampler = Metropolis(propose)
  d.precisionInvScale = sample(sampler, logratio, d.precisionInvScale)
end
