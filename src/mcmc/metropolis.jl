import StatsBase.sample

export
    Metropolis,
    sample

mutable struct Metropolis
  propose::Function
  acceptancerate::Float64
  num::Float64
end

Metropolis(propose::Function) = Metropolis(propose, 0.0, 0.0)

function sample(sampler::Metropolis, logratio::Function, value::Float64)
  newvalue = sampler.propose(value)
  arate = min(1.0, exp(logratio(newvalue,value)))
  sampler.num += 1
  sampler.acceptancerate += (arate - sampler.acceptancerate) / sampler.num
  if rand() < arate
    return newvalue
  else
    return value
  end
end
