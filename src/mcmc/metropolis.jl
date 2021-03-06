export
    Metropolis,
    sample

mutable struct Metropolis
  propose         ::  Function  # Proposal function
  acceptancerate  ::  Float     # Acceptance rate
  num             ::  Float     # Number of iterations
end

Metropolis(propose::Function) = Metropolis(propose, 0.0, 0.0)

function sample(sampler::Metropolis, logratio::Function, value::Float)
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
