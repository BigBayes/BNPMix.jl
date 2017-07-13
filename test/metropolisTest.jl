include("../src/mcmc/metropolis.jl")
using Base.Test, Distributions

numdata = 500000

precisionInvScaleAlpha = 0.2
precisionInvScaleBeta = 6.34
precisionInvScale = precisionInvScaleAlpha / precisionInvScaleBeta

function propose(x::Float64)
  x*exp(.5*rand(Normal()))
end

function logratio(isnew::Float64, isold::Float64)
  return precisionInvScaleAlpha * (log(isnew) - log(isold)) - precisionInvScaleBeta * (isnew - isold)
end

sampler = Metropolis(propose)

mean = 0.0
var = 0.0

for i in 1:numdata
  precisionInvScale = sample(sampler, logratio, precisionInvScale)
  mean += precisionInvScale
  var += precisionInvScale*precisionInvScale
end

mean /= numdata
var = var / (numdata - 1) - mean*mean

trueMean = 0.03
trueVar = 0.005

ε = 0.01
@test ≈(mean, trueMean, atol=ε)
@test ≈(var, trueVar, atol=ε)
