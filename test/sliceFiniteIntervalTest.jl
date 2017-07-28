using BNPMix
using Base.Test

sampler = SliceFiniteInterval(0.0, 1.0)
sigma = 0.5
numdata = 10000;

function logDensity(x::Float64) #Beta(2,5)
 return log(x) + 4*log(1-x)
end

mean = 0.0
var = 0.0

for i in 1:numdata
  sigma = sample(sampler, logDensity, sigma)
  mean += sigma
  var += sigma*sigma
end

mean /= numdata
var = var / (numdata - 1) - mean*mean

trueMean = 2.0 / 7.0
trueVar = 10.0 / (49.0 * 8.0)

ε = 0.01
@test ≈(mean, trueMean, atol=ε)
@test ≈(var, trueVar, atol=ε)
