using BNPMix
using Base.Test

sampler = SliceStepOut(0.2, 50)
sigma = 0.5
numdata = 10000

function logDensity(x::Float64)
 return - (x+1)^2/(2*4) - (x-1)^2/(2*4)
end

mean = 0.0
var = 0.0

for i in 1:numdata
  #sample(sampler, logDensity, sigma)
  sigma = sample(sampler, logDensity, sigma)
  mean += sigma
  var += sigma*sigma
end

mean /= numdata
var = var / (numdata - 1) - mean*mean

trueMean = 0.0
trueVar = 2.0

ε = 0.1
@test ≈(mean, trueMean, atol=ε)
@test ≈(var, trueVar, atol=ε)
