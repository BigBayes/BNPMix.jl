using Distributions
import StatsBase.sample
include("../mcmc/sliceStepOut.jl")
include("../mixture/cluster.jl")

# Normalized Gamma Generalized Process
abstract type NRMI end

mutable struct NGGP <: NRMI
  Ashape::Float64
  Ainvscale::Float64
  Salpha::Float64
  Sbeta::Float64
  Tshape::Float64
  Tinvscale::Float64
  alpha::Float64
  sigma::Float64
  tau::Float64
  logU::Float64
end

NGGP(Ashape::Float64, Ainvscale::Float64, Salpha::Float64, Sbeta::Float64, Tshape::Float64, Tinvscale::Float64) =
NGGP(Ashape, Ainvscale, Salpha, Sbeta, Tshape, Tinvscale, (Ashape+1.0)/(Ainvscale+1.0), 0.1, (Tshape+1.0)/(Tinvscale+1.0), 0.0)

function logLevy{T<:NRMI}(mu::T, mass::Float64)
  logLevy(mu, mass, 0.0)
end

function logLevy(nggp::NGGP, mass::Float64, logu::Float64)
  logLevy(mass, nggp.alpha, nggp.sigma, nggp.tau, logu)
end

function logLevy(mass::Float64, alpha::Float64, sigma::Float64, tau::Float64, logu::Float64)
  return log(alpha) - lgamma(1.0-sigma) - (1.0+sigma)*log(mass) - (exp(logu)+tau)*mass
end

function laplace(nggp::NGGP)
  return laplace(nggp.alpha, nggp.sigma, nggp.tau, nggp.logU)
end

function laplace(nggp::NGGP, logu::Float64)
  return laplace(nggp.alpha, nggp.sigma, nggp.tau, logu)
end

function laplace(alpha::Float64, sigma::Float64, tau::Float64, logu::Float64)
  if (sigma<1e-16) return alpha*log1p(exp(logu)/tau)
	else return alpha/sigma*((tau+exp(logu))^sigma-tau^sigma)
  end
end

function logGamma(nggp::NGGP, num::Int64)
  return logGamma(nggp.alpha, nggp.sigma, nggp.tau, num, nggp.logU)
end

function logGamma(nggp::NGGP, num::Int64, logu::Float64)
  return logGamma(nggp.alpha, nggp.sigma, nggp.tau, num, logu)
end

function logGamma(alpha::Float64, sigma::Float64, tau::Float64, num::Int64, logu::Float64)
  return log(alpha) - lgamma(1.0-sigma) + lgamma(num-sigma) - (num-sigma)*log(tau+exp(logu))
end

function logMeanMass(nggp::NGGP, num::Int64)
  return  logMeanMass(nggp.sigma, nggp.tau, num, nggp.logU)
end

#function logMeanMass(nggp::NGGP, num::Int64)
  #return  logMeanMass(nggp, num, nggp.logU)
#end

function logMeanMass(sigma::Float64, tau::Float64, num::Int64, logu::Float64)
  return log((num-sigma) / (tau+exp(logu)))
end

#function logMeanMass(nggp::NGGP, num::Int64, logu::Float64)
  #return log(num-nggp.sigma) - log(nggp.tau+exp(logu))
  #return log((num-nggp.sigma) / (nggp.tau+exp(logu)))
#end

function logMeanTotalMass(nggp::NGGP, numclusters::Int64) # why numclusters ??
  return logMeanTotalMass(nggp.alpha, nggp.sigma, nggp.tau, numclusters, nggp.logU)
end

function logMeanTotalMass(alpha::Float64, sigma::Float64, tau::Float64, numclusters::Int64, logu::Float64) # why numclusters ??
  return log(alpha) + (sigma-1.0)*log(tau+exp(logu))
end

function drawLogMass{T<:NRMI}(mu::T, num::Int64)
  return drawLogMass(mu, num, mu.logU)
end

function drawLogMass(nggp::NGGP, num::Int64, logu::Float64)
  return log(rand(Gamma(num-nggp.sigma,1))) - log(nggp.tau+exp(logu))
end

function meanNumClusters(nggp::NGGP, num::Int64)
  return meanNumClusters(nggp, num, nggp.logU)
end

function meanNumClusters(nggp::NGGP, num::Int64, logu::Float64)
  print("NGGP.meanNumClusters: dont know answer yet\n")
  return nggp.alpha*log(1+num/nggp.alpha)
end

### Sampling
function sample(nggp::NGGP, numData::Int64, clusters::Set{Cluster})
  #print("NGGP - sample\n")
  sampleU(nggp, numData, clusters)
  sampleAlpha(nggp, clusters)
  sampleSigma(nggp, clusters)
  sampleTau(nggp, clusters)
end

function sampleU(nggp::NGGP, numData::Int64, clusters::Set{Cluster}) # No prior on hyperparameters
  #print("NGGP - sampleU\n")
  sampler = SliceStepOut(1.0, 20)

  function logDensityU(logu::Float64)
    result = numData *  logu - laplace(nggp, logu)
    for cc in clusters
      assert(!isEmpty(cc))
      result += logGamma(nggp, cc.number, logu)
    end
    return result
  end

  nggp.logU = sample(sampler, logDensityU, nggp.logU)
end

function sampleAlpha(nggp::NGGP, clusters::Set{Cluster})
  shape = nggp.Ashape + length(clusters)
  invscale = nggp.Ainvscale +
    ((nggp.sigma<1e-16) ? log1p(exp(nggp.logU)/nggp.tau) :
    ((nggp.tau+exp(nggp.logU))^nggp.sigma-(nggp.tau)^nggp.sigma)/nggp.sigma)
  gamma = Gamma(shape, 1.0/invscale)
  nggp.alpha = rand(gamma)
end

function sampleSigma(nggp::NGGP, clusters::Set{Cluster})
  sampler = SliceFiniteInterval(0.0, 1.0)

  function logDensity(x::Float64)
    result = (nggp.Salpha-1.0)*log(max(1e-16,x)) + (nggp.Sbeta-1.0)*log(max(1e-16,1.0-x)) -
       laplace(nggp.alpha, x, nggp.tau, nggp.logU)
    for cc in clusters
      if !isEmpty(cc)
        result += logGamma(nggp.alpha, x, nggp.tau, cc.number, nggp.logU)
      end
    end
    return result
  end

  nggp.sigma = sample(sampler, logDensity, nggp.sigma)
end

function sampleTau(nggp::NGGP, clusters::Set{Cluster})
  function propose(x::Float64)
    x*exp(.5*rand(Normal()))
  end

  function logratio(Tnew::Float64, Told::Float64)
    result = nggp.Tshape*log(Tnew/Told) - nggp.Tinvscale*(Tnew-Told) -
      laplace(nggp.alpha, nggp.sigma, Tnew, nggp.logU) +
      laplace(nggp.alpha, nggp.sigma, Told, nggp.logU)

      for cc in clusters
        if !isEmpty(cc)
          result = result + logGamma(nggp.alpha, nggp.sigma, Tnew, cc.number, nggp.logU) -
                            logGamma(nggp.alpha, nggp.sigma, Told, cc.number, nggp.logU)
        end
      end
      return result
  end

  sampler = Metropolis(propose)
  nggp.tau = sample(sampler, logratio, nggp.tau)
end
