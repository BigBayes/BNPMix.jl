import Base.run
using ProgressMeter

export
    Sampler,
    runSampler

struct Sampler
  model::Mixture
  numBurnIn::Int
  numSample::Int
  numThinning::Int
  outputFilename::String
  collectors::Array{Float}{2}
end

Sampler(model::Mixture, numBurnIn::Int, numSample::Int, numThinning::Int, outputFilename::String) =
Sampler(model, numBurnIn, numSample, numThinning, outputFilename, zeros(Float, numSample, 6))

function runSampler(s::Sampler)
  print("\n Initialize Sampler: \n")
  initializeSampler(s.model)
  p = Progress(s.numBurnIn, 1, "Burn-in: ", 50)
  tic()
  for i in 1:s.numBurnIn
    sample(s.model)
    update!(p, i)
  end
  lruntime = toq()
  #Flush
  p = Progress(s.numSample, 5, "Sampling: ", 50)
  for i in 1:s.numSample
    tic()
    sample(s, s.numThinning)
    lruntime += toq()
    update!(p, i)
    s.collectors[i, 1] = length(s.model.clusters)
    s.collectors[i, 2] = s.model.nrmi.alpha
    s.collectors[i, 3] = s.model.nrmi.sigma
    s.collectors[i, 4] = s.model.nrmi.tau
    s.collectors[i, 5] = s.model.nrmi.logU
    s.collectors[i, 6] = s.model.prior.precisionInvScale
    #Collect
    #Flush
  end
  print("\n Done: ", "\nRun time = ", lruntime)
  outputFilename = s.outputFilename
  Base.run(`rm -f output/$outputFilename`)
  Base.run(`touch output/$outputFilename`)
  open(string("output/", outputFilename), "w") do f
     for i in 1:s.numSample
        n1 = s.collectors[i, 1]
        n2 = s.collectors[i, 2]
        n3 = s.collectors[i, 3]
        n4 = s.collectors[i, 4]
        n5 = s.collectors[i, 5]
        n6 = s.collectors[i, 6]
        write(f, "$n1 $n2 $n3 $n4 $n5 $n6\n")
     end
   end
end

function sample(s::Sampler, numIteration::Int)
  for i in 1:numIteration
    sample(s.model)
  end
end
