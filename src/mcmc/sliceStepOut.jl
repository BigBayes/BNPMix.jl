include("../utils.jl")
using Distributions
import Base.convert

export
    SliceFiniteInterval,
    SliceStepOut,
    sample

### SliceFiniteInterval

struct SliceFiniteInterval
  lower::Float
  upper::Float
end

function sample(sampler::SliceFiniteInterval, logDensity::Function, value::Float)
   slice = logDensity(value) - rndExponential()
   return  sample(sampler, logDensity, value, slice)
end


function sample(sampler::SliceFiniteInterval, logDensity::Function, value::Float, slice::Float)
  #print("SliceFiniteInterval - sample\n")
  l = sampler.lower
  u = sampler.upper
  while true
    newvalue = rand(Uniform(l, u))
    if logDensity(newvalue)>slice return newvalue
    elseif newvalue>value u = newvalue
    else l = newvalue
    end
  end
end

### SliceStepOut

mutable struct SliceStepOut
  stepsize::Float
  numstep::Int
  lower::Float
  upper::Float
end

SliceStepOut(stepsize::Float, numstep::Int) = SliceStepOut(stepsize, numstep, 0.0, 1.0)

convert(::Type{SliceFiniteInterval}, sampler::SliceStepOut) = SliceFiniteInterval(sampler.lower, sampler.upper)

function setLower(sampler::Union{SliceFiniteInterval, SliceStepOut}, l::Float)
  sampler.lower = l
end

function setUpper(sampler::Union{SliceFiniteInterval, SliceStepOut}, u::Float)
  sampler.upper = u
end

function sample(sampler::SliceStepOut, logDensity::Function, value::Float)
  #print("SliceStepOut - sample\n")
  slice = logDensity(value) - rndExponential()
  l = value - sampler.stepsize * rand(Uniform())
  u = l + sampler.stepsize
  nl = round(rand(Uniform())*sampler.numstep)
  for i in 0:nl-1
    y = logDensity(l)
    if y < slice
      break
    end
    l -= sampler.stepsize
  end
  nu = sampler.numstep - 1 - nl
  for i in 0:nu-1
    y = logDensity(u)
    if y < slice
      break
    end
    u += sampler.stepsize
  end
  setLower(sampler, l)
  setUpper(sampler, u)

  return sample(convert(SliceFiniteInterval, sampler), logDensity, value, slice)
end
