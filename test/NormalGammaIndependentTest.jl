using BNPMix
using Base.Test, Distributions

meanMean = 2.17255
meanPrecision = 0.634557
precisionShape = 2.0
prior = NormalGammaIndependent(meanMean, meanPrecision, precisionShape, 0.2, 6.34557)

numdata = 100000
mean = 0.0
var = 0.0

for i in 1:numdata
  normal = drawSample(prior)
  mean += normal.μ
  var += normal.σ
end

mean /= numdata
var /= numdata

ε = 0.05
@test ≈(mean, meanMean, atol=ε)
@test ≈(var, 1/sqrt(precisionShape/prior.precisionInvScale), atol=ε)


@test ≈(logProbability(prior, Normal(1.0, 1/sqrt(0.1))), -10.81, atol=ε)
