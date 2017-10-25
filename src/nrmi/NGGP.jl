export
    NGGP,
    logLevy,
    laplace,
    logGamma,
    logMeanMass,
    logMeanTotalMass,
    meanNumClusters

abstract type NRMI end

mutable struct NGGP <: NRMI # Normalized Gamma Generalized Process
  Ashape     ::  Float  # Gamma's shape parameter for alpha prior
  Ainvscale  ::  Float  # Gamma's inverse scale parameter for alpha prior
  Salpha     ::  Float  # Beta's alpha parameter for sigma prior
  Sbeta      ::  Float  # Beta's beta parameter for sigma prior
  Tshape     ::  Float  # Gamma's shape parameter for tau prior
  Tinvscale  ::  Float  # Gamma's inverse scale parameter for tau prior
  alpha      ::  Float  # alpha parameter of NGGP's Lévy intensity measure
  sigma      ::  Float  # sigma parameter of NGGP's Lévy intensity measure
  tau        ::  Float  # tau parameter of NGGP's Lévy intensity measure
  logU       ::  Float  # Auxiliary random variable
end

NGGP(Ashape::Float, Ainvscale::Float, Salpha::Float, Sbeta::Float, Tshape::Float, Tinvscale::Float) =
NGGP(Ashape, Ainvscale, Salpha, Sbeta, Tshape, Tinvscale, (Ashape+1.0)/(Ainvscale+1.0), 0.1, (Tshape+1.0)/(Tinvscale+1.0), 0.0)

NGGP(alpha::Float, sigma::Float, tau::Float) = new(.0,.0,.0,.0,.0,.0,alpha,sigma,tau,0.0)

function logLevy{T<:NRMI}(mu::T, mass::Float)
  logLevy(mu, mass, 0.0) #logLevy(mu, mass, mu.logU) ??
end

function logLevy(nggp::NGGP, mass::Float, logu::Float)
  logLevy(mass, nggp.alpha, nggp.sigma, nggp.tau, logu)
end

function logLevy(mass::Float, alpha::Float, sigma::Float, tau::Float, logu::Float)
  return log(alpha) - lgamma(1.0-sigma) - (1.0+sigma)*log(mass) - (exp(logu)+tau)*mass
end

function laplace(nggp::NGGP)
  return laplace(nggp.alpha, nggp.sigma, nggp.tau, nggp.logU)
end

function laplace(nggp::NGGP, logu::Float)
  return laplace(nggp.alpha, nggp.sigma, nggp.tau, logu)
end

function laplace(alpha::Float, sigma::Float, tau::Float, logu::Float)
  if (sigma<1e-16) return alpha*log1p(exp(logu)/tau)
	else return alpha/sigma*((tau+exp(logu))^sigma-tau^sigma)
  end
end

function logGamma(nggp::NGGP, num::Int)
  return logGamma(nggp.alpha, nggp.sigma, nggp.tau, num, nggp.logU)
end

function logGamma(nggp::NGGP, num::Int, logu::Float)
  return logGamma(nggp.alpha, nggp.sigma, nggp.tau, num, logu)
end

function logGamma(alpha::Float, sigma::Float, tau::Float, num::Int, logu::Float)
  return log(alpha) - lgamma(1.0-sigma) + lgamma(num-sigma) - (num-sigma)*log(tau+exp(logu))
end

function logMeanMass(nggp::NGGP, num::Int)
  return  logMeanMass(nggp.sigma, nggp.tau, num, nggp.logU)
end

function logMeanMass(sigma::Float, tau::Float, num::Int, logu::Float)
  return log((num-sigma) / (tau+exp(logu)))
end

function logMeanTotalMass(nggp::NGGP, numclusters::Int) # why numclusters ??
  return logMeanTotalMass(nggp.alpha, nggp.sigma, nggp.tau, numclusters, nggp.logU)
end

function logMeanTotalMass(alpha::Float, sigma::Float, tau::Float, numclusters::Int, logu::Float) # why numclusters ??
  return log(alpha) + (sigma-1.0)*log(tau+exp(logu))
end

function drawLogMass{T<:NRMI}(mu::T, num::Int)
  return drawLogMass(mu, num, mu.logU)
end

function drawLogMass(nggp::NGGP, num::Int, logu::Float)
  return log(rand(Gamma(num-nggp.sigma,1))) - log(nggp.tau+exp(logu))
end

function meanNumClusters(nggp::NGGP, num::Int)
  return meanNumClusters(nggp, num, nggp.logU)
end

function meanNumClusters(nggp::NGGP, num::Int, logu::Float)
  return nggp.alpha*log(1+num/nggp.alpha)
end

### Sampling
function sample(nggp::NGGP, numData::Int, clusters::Set{Cluster})
  sampleU(nggp, numData, clusters)
  sampleAlpha(nggp, clusters)
  sampleSigma(nggp, clusters)
  sampleTau(nggp, clusters)
end

function sampleU(nggp::NGGP, numData::Int, clusters::Set{Cluster}) # No prior on hyperparameters
  sampler = SliceStepOut(1.0, 20)

  function logDensityU(logu::Float)
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

  function logDensity(x::Float)
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
  function propose(x::Float)
    x*exp(.5*rand(Normal()))
  end

  function logratio(Tnew::Float, Told::Float)
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

function logW(nggp::NGGP, s::Float, t::Float)
  return log(nggp.alpha) -
  lgamma(1.0-nggp.sigma) - (1.0+nggp.sigma)*log(t) - s*(exp(nggp.logU)+nggp.tau)
end

function invWInt(nggp::NGGP, x::Float, t::Float, logu::Float)
  a = log(x*(exp(logu)+nggp.tau))
  b = logLevy(nggp, t, logu)
  if a >= b  return Inf end
  if b > a+33.0 return t+exp(a-b)/(exp(logu)+nggp.tau)
  else return t-log(1.0-exp(a-b))/(exp(logu)+nggp.tau)
  end
end

function drawLogMasses(nggp::NGGP, slice::Float)
  minSlice = 1.0e-7
  maxClusters = 1000000
  slice = exp(slice)
  if slice < minSlice slice = 1.0e-6 end

  s = slice
  masses = Array{Float}(0)
  while(true)
    e = rand(Exponential(1.0))
    news = invWInt(nggp, e, s, nggp.logU)
    if news == Inf break end
    if (rand() < exp(logLevy(nggp, news, nggp.logU) - logW(nggp, news, s)))
      push!(masses, news)
    end
    if length(masses) > maxClusters
      masses = Array{Float}(0)
      slice *= 10.0
      news = slice
    end
    s = news
  end

  masses = masses[end:-1:1]
  result = zeros(Float, length(masses))
  for i in 1:length(masses)
    result[i] = log(masses[i])
  end

  return result
end
