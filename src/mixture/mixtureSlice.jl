export
    MixtureSlice

mutable struct MixtureSlice <: Mixture
  nrmi              ::  NRMI    # Normalized random measure
  prior             ::  Prior   # Base distribution (prior over component parameters)
  factory           ::  Factory # Factory object to generate component parameters from prior
  data              ::  Union{Void, Array{Float, 2}}  # Observed data
  map               ::  Union{Array{Union{Cluster, Void}, 0}, Array{Union{Cluster, Void}, 1}} # Mapping from data to clusters
  clusters          ::  Set{Cluster} # Clusters making up the mixture
  numClustersRatio  ::  Float
  minslice          ::  Float
  slice             ::  Union{Array{Float, 0},Array{Float, 1}}
  clusterarray      ::  Union{Array{Cluster, 0},Array{Cluster, 1}}

  MixtureSlice(nrmi::NRMI, prior::Prior, factory::Factory) =
    new(nrmi, prior, factory, nothing, Array{Union{Void, Cluster}}(0), Set{Cluster}(), 5.0, 1.0e-8, Array{Float}(0), Array{Cluster}(0))
end

function clearClusterarray(m::MixtureSlice)
  m.clusterarray = Array{Cluster}(0)
end

function addClustersToClusterarray(m::MixtureSlice)
  for cc in m.clusters
    push!(m.clusterarray, cc)
  end
end

function addData(m::MixtureSlice, traindata::Array{Float, 2})
  addDataAbstract(m, traindata)

  numdata = size(traindata, 1)
  for i in 1:numdata
    push!(m.slice, getCluster(m, i).logmass - rndExponential())
  end
  # update partition structure variables, but not cluster parameters
  for cc in m.clusters
    cc.logmass = drawLogMass(m.nrmi, cc.number)
    #sample(cc.parameter)
  end
  minMass = Inf
  for i in 1:numdata
    s = getCluster(m, i).logmass - rndExponential()
    m.slice[i] = s
    minMass = min(minMass, s)
  end

  masses = drawLogMasses(m.nrmi, minMass)
  clearClusterarray(m)
  addClustersToClusterarray(m)
  for mass in masses
    push!(m.clusterarray, newCluster(m, mass))
  end
  m.clusterarray = sort(m.clusterarray, by=cc->cc.logmass, rev=true)

  sampleAssignments(m)
end

function sampleClusters(m::MixtureSlice)
  for cc in m.clusters
    cc.logmass = drawLogMass(m.nrmi, cc.number)
    sample(cc.parameter)
  end

  minMass = Inf
  for i in 1:size(m.data,1)
    s = getCluster(m, i).logmass - rndExponential()
    m.slice[i] = s
    minMass = min(minMass, s)
  end

  masses = drawLogMasses(m.nrmi, minMass)
  clearClusterarray(m)
  addClustersToClusterarray(m)
  for mass in masses
    push!(m.clusterarray, newCluster(m, mass))
  end
  m.clusterarray = sort(m.clusterarray, by=cc->cc.logmass, rev=true)
end

function sampleAssignments(m::MixtureSlice)
  numData = size(m.data, 1)
  for i in 1:numData
    sampleAssignment(m, i)
  end
  m.clusters = Set{Cluster}()
  num = 0
  for cc in m.clusterarray
    if !isEmpty(cc)
      push!(m.clusters, cc)
    end
    num += cc.number
  end
  assert(num == numData)
end

function sampleAssignment(m::MixtureSlice, index::Int)
  datum = getDatum(m, index)
  unassign(m, index, datum)
  thisslice = m.slice[index]

  for cc in m.clusterarray
    if cc.logmass < thisslice break end
    cc.w = logPredictive(cc.parameter, datum)
  end

  ma = -Inf
  for cc in m.clusterarray
    if cc.logmass < thisslice break end
    if ma < cc.w ma = cc.w end
  end

  s = 0.0
  for cc in m.clusterarray
    if cc.logmass < thisslice break end
    cc.w = exp(cc.w-ma)
    s += cc.w
  end

  r = rand(Uniform())*s
  for cc in m.clusterarray
    r -= cc.w
    if r <= 0.0
      assign(m, index, datum, cc)
      break
    end
  end

  assert(r <= 0.0)
end
