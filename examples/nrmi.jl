using BNPMix

# alg parameters
numEmptyClusters = 2
numBurnin = 10000
numSample = 10000
numThinning = 20
conjugate = false
alg = "Neal8"#"Reuse"

function nrmi(alg::String, conjugate::Bool, outputFilename::String, numBurnin::Int64, numSample::Int64, numThinning::Int64, numEmptyClusters::Int64)

  numdim = size(data, 2)
  numdata = size(data, 1)

  # distributions hyperparameters
  alphaShape = 1.0
  alphaInvScale = 1.0
  sigmaAlpha = 1.0
  sigmaBeta = 2.0
  tauShape = 1e9
  tauInvScale = 1e9

  meanRelScale = 1.0
  precisionScale = 50.0
  invScaleDegFreedom = numdim-0.6
  precisionDegFreedom = numdim+3.0
  invScaleInvScale = precisionScale*invScaleDegFreedom/(precisionDegFreedom-numdim-1.0)

  # Compute elementary statistics

  dmin = ones(Float64,numdim) * Inf
  dmax = - ones(Float64,numdim) * Inf
  dmean = zeros(Float64,numdim)
  dvar = zeros(Float64,numdim)

  for d in 1:numdim
    for i in 1:numdata
      x = data[i][d]
      dmean[d] += x
      dvar[d] += x*x
      if x < dmin[d] dmin[d] = x
      elseif x > dmax[d] dmax[d] = x
      end
    end
    dmean[d] /= numdata
    dvar[d] = dvar[d] / (numdata-1) - dmean[d] * dmean[d]
  end

  meanMean = zeros(Float64, numdim)
  meanPrecision = zeros(Float64, numdim, numdim)
  precisionInvScaleInvScale = zeros(Float64, numdim, numdim)

  for d in 1:numdim
    meanMean[d] = 0.5 * (dmax[d] + dmin[d])
    range = 0.5 * (dmax[d] - dmin[d])
    meanPrecision[d,d] = 1.0 / range / range / meanRelScale
    precisionInvScaleInvScale[d,d] = invScaleInvScale / range / range
  end

  # initialize prior and factory

  if numdim > 1
    print("Multidimensional not implemented\n")
    #MVNormalWishartIndependent
  else
    if conjugate
      print("Conjugate not implemented\n")
    else
      meanmean = meanMean[1]
      meanprecision = meanPrecision[1,1]
      precisioninvscaleinvscale = precisionInvScaleInvScale[1,1]
      prior = NormalGammaIndependent(meanmean, meanprecision, .5*precisionDegFreedom, .5*invScaleDegFreedom, precisioninvscaleinvscale)
      factory = NormalNonConjugateFactorySampled()
    end
  end

  # build model

  nggp = NGGP(alphaShape, alphaInvScale ,sigmaAlpha, sigmaBeta, tauShape, tauInvScale)

  if alg == "Neal8" model = MixtureNeal8(nggp, prior, factory, numEmptyClusters)
  elseif alg == "Reuse" model = MixtureReuse(nggp, prior, factory, numEmptyClusters)
  elseif alg == "Slice" model = MixtureSlice(nggp, prior, factory)
  else print("Only Neal8, Reuse and Slice implemented\n")
  end

  if numdim == 1
    newData = zeros(numdata, 1)
    newData[:,1] = data
  else
    newData = data
  end

  addData(model, newData)

  # Parameters STDOUT
  print("\nNRMIX ", alg);
  print("\n  components: ", conjugate ? "conjugate" : "non-conjugate")
  print("\n  outputFilename = ", outputFilename)
  print("\n  numDimension = ", numdim)
  print("\n  numData = ", numdata)
  print("\n  numBurnin = ", numBurnin)
  print("\n  numSample = ", numSample)
  print("\n  numThinning = ", numThinning)
  print("\n  numNewClusters = ", numEmptyClusters)
  print("\n  alphaShape = ", alphaShape)
  print("\n  alphaInvScale = ", alphaInvScale)
  print("\n  sigmaAlpha = ", sigmaAlpha)
  print("\n  sigmaBeta = ", sigmaBeta)
  print("\n  tauShape = ", tauShape)
  print("\n  tauInvScale = ", tauInvScale)
  print("\n  meanRelativeScale = ", meanRelScale)
  print("\n  precisionDegFreedom = ", precisionDegFreedom)
  print("\n  invScaleDegFreedom = ", invScaleDegFreedom)
  print("\n  precisionScale = ", precisionScale)
  #print("\n  useMeanVar = ", useMeanVar)

  # Sampling
  collectors = ["nbClusters", "nrmi_sigma", "nrmi_tau", "nrmi_logU", "prior_precisionInvScale"]
  sampler = Sampler(model, collectors, numBurnin, numSample, numThinning, outputFilename)
  runSampler(sampler)

end
