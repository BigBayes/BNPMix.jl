export
    MixtureReuse,
    addData

mutable struct MixtureReuse <: Mixture
  nrmi              ::  NRMI    # Normalized random measure
  prior             ::  Prior   # Base distribution (prior over component parameters)
  factory           ::  Factory # Factory object to generate component parameters from prior
  numNewClusters    ::  Int     # Number of clusters making up the augmented space
  data              ::  Union{Void, Array{Float, 2}}  # Observed data
  map               ::  Union{Array{Union{Cluster, Void}, 0}, Array{Union{Cluster, Void}, 1}} # Mapping from data to clusters
  clusters          ::  Set{Cluster} # Clusters making up the mixture
  numClustersRatio  ::  Float
  empties           ::  Array{Union{Cluster, Void}, 1} # Augmented space

  function MixtureReuse(nrmi::NRMI, prior::Prior, factory::Factory, numNewClusters::Int)
    this = new(nrmi, prior, factory, numNewClusters, nothing, Array{Union{Void, Cluster}}(0), Set{Cluster}(), 5.0, Array{Union{Void, Cluster}}(0))
    initializeEmpties(this)
    return this
  end
end

function initializeEmpties(m::MixtureReuse)
  for i in 1:m.numNewClusters
    cc = newCluster(m, 0.0)
    push!(m.empties, cc)
  end
end

function addData(m::MixtureReuse, traindata::Array{Float, 2})
  addDataAbstract(m, traindata)
  sampleAssignments(m)
end

function sampleAssignments(m::MixtureReuse)
  numData = size(m.data, 1)
  for cc in m.clusters
    cc.logmass = logMeanMass(m.nrmi, cc.number)
  end
  for i in 1:m.numNewClusters
    cc = m.empties[i]
    cc.logmass = nothing
    sample(cc.parameter)
    push!(m.clusters, cc)
  end
  for i in 1:numData
    sampleAssignment(m, i)
  end
  for i in 1:m.numNewClusters
    cc = m.empties[i]
    assert(isEmpty(cc))
    delete!(m.clusters, cc)
  end
end

function sampleAssignment(m::MixtureReuse, index::Int)
  datum = getDatum(m, index)
  tt = unassign(m, index, datum)

  # initialize empty clusters
  if !isEmpty(tt)
    tt.logmass = logMeanMass(m.nrmi, tt.number)
  else
    k = rand(1:m.numNewClusters)
    newk = m.empties[k]
    #factory.destruct(prior, newk.parameter) # useless
    newk.parameter = tt.parameter
    delete!(m.clusters, tt)
  end

  # compute probabilities
  mi = -Inf
  ew = logMeanTotalMass(m.nrmi, length(m.clusters)-m.numNewClusters) - log(m.numNewClusters)
  for cc in m.clusters
    cm = cc.logmass
    if (cm == nothing) cc.w = ew
    else cc.w = cm
    end
    cc.w += logPredictive(cc.parameter, datum)
    if mi < cc.w mi = cc.w
    end
  end

  s = 0.0
  for cc in m.clusters
    cc.w = exp(cc.w-mi)
    s += cc.w
  end

  # sample assignment
  r = rand(Uniform())*s
  for cc in m.clusters
    r -= cc.w
    if r <= 0.0
      if !isEmpty(cc)
        assign(m, index, datum, cc)
        cc.logmass = logMeanMass(m.nrmi, cc.number)
      else
        # create new cluster and copy parameter over
        tt = newCluster(m, logMeanMass(m.nrmi, 1))
        param = tt.parameter
        tt.parameter = cc.parameter
        cc.parameter = param
        push!(m.clusters, tt)
        assign(m, index, datum, tt)
      end
      break
    end
  end

  assert(r <= 0.0)
end
