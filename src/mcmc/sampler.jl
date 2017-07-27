export
    Sampler,
    runSampler

struct Sampler
  model           ::  Mixture
  collectors      ::  Array{Any}{1}  # Array of variables to be collected
  numBurnIn       ::  Int            # Number of burn-in iterations
  numSample       ::  Int            # Number of samples
  numThinning     ::  Int            # Number of thinning iterations
  outputFilename  ::  String
end

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

  p = Progress(s.numSample, 5, "Sampling: ", 50)
  collector = zeros(Float, s.numSample, length(s.collectors))
  for i in 1:s.numSample
    tic()
    sample(s, s.numThinning)
    lruntime += toq()
    update!(p, i)

    for j in 1:length(s.collectors)
      collector[i, j] = get(s.model, s.collectors[j])
    end
  end

  print("\n Done: ", "\nRun time = ", lruntime)
  outputFilename = s.outputFilename
  Base.run(`mkdir -p output`)
  Base.run(`rm -f output/$outputFilename`)
  Base.run(`touch output/$outputFilename`)
  open(string("output/", outputFilename), "w") do f
     for i in 1:s.numSample
       line = string()
       for j in 1:length(s.collectors)
         line = string(line, " ",  collector[i, j])
       end
        write(f, "$line\n")
     end
   end
end

function sample(s::Sampler, numIteration::Int)
  for i in 1:numIteration
    sample(s.model)
  end
end
