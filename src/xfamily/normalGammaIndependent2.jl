using Distributions

export
    NormalGammaIndependent,
    drawSample,
    logProbability

abstract type Hierarchy end
abstract type Prior end

mutable struct NormalGammaIndependent <: Prior
    meanMean                ::  Float  # Hyperparameters of Normal-Gamma-Independent
    meanPrecision           ::  Float  # idem
    precision               ::  Float  # idem
    data                    ::  Union{Set{Hierarchy}, Void}

    function NormalGammaIndependent(meanMean::Float, meanPrecision::Float, precision::Float)
    	new(meanMean, meanPrecision, precision, nothing)
    end
end

function drawSample(d::NormalGammaIndependent)
  return Normal(rand(Normal(d.meanMean, 1.0/sqrt(d.meanPrecision))), 1.0/sqrt(d.precision))
end

function logNormalizer(d::NormalGammaIndependent)
  halfLogTwoPi = .5*log(2*pi)
  return halfLogTwoPi -.5*log(d.meanPrecision)
end

function logProbability(d::NormalGammaIndependent, normal::Distributions.Normal{Float})
  diff = mean(normal) - d.meanMean
  normalprecision = 1.0/var(normal)
	return  -.5*d.meanPrecision*diff*diff +
					logNormalizer(d)
end

### Sampling

function sample(d::NormalGammaIndependent, data::Set{Hierarchy})
  d.data = data
  #sampleMeanMean(d)
  #sampleMeanPrecision(d)
end
