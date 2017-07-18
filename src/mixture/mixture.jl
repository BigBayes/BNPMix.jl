abstract type Mixture end

export
    assign,
    unassign,
    getDatum

# Helpers
function numData{T<:Mixture}(m::T)
  return size(m.data, 1)
end

function numClusters{T<:Mixture}(m::T)
  return length(m.clusters)
end

function getDatum{T<:Mixture}(m::T, i::Int)
  assert((i >= 1) & (i <= length(m.map)))
  return m.data[i,:]
end

function getCluster{T<:Mixture}(m::T, i::Int)
  assert((i >= 1) & (i <= length(m.map)))
  return m.map[i]
end

function newCluster{T<:Mixture}(m::T, logmass::Float)
  return Cluster(0, logmass, construct(m.factory, m.prior))
end

# Clusters management
function addDataAbstract{T<:Mixture}(m::T, traindata::Array{Float, 2})
  m.data = traindata
  numdata = numData(m)
  numclusters = Int(min(numdata,round(m.numClustersRatio*ceil(meanNumClusters(m.nrmi,numdata)))))
  clusterlist = Array{Cluster}(0)
  @inbounds for i in 1:numclusters
    cc = newCluster(m, logMeanMass(m.nrmi, 1))
    push!(clusterlist, cc)
    push!(m.clusters, cc)
    push!(m.map,nothing)
    assign(m, i, m.data[i,:], cc)
  end
  @inbounds for i in numclusters+1:numdata
    cc = clusterlist[rand(1:numclusters)]
    push!(m.map,nothing)
    assign(m, i, m.data[i,:], cc)
  end
end

function assign{T<:Mixture}(m::T, i::Int, datum::Array{Float, 1}, cc::Cluster)
  assert((length(m.map)>=i) | (m.map[i] == nothing))
  assert(datum == getDatum(m, i))
  addDatum(cc.parameter, datum)
  cc.number += 1
  m.map[i] = cc
end

function unassign{T<:Mixture}(m::T, i::Int, datum::Array{Float, 1})
  assert(((datum == nothing) & (getDatum(m, i)==nothing)) | (datum == getDatum(m, i)))
  assert((i>=1) & (i<=length(m.map)))
  cc = m.map[i]
  m.map[i] = nothing
  assert(cc != nothing)
  assert(cc.parameter != nothing)
  if datum != nothing
    removeDatum(cc.parameter, datum)
  end
  cc.number -= 1
  assert(cc.number >= 0)
  return cc
end

# Sampling

function initializeSampler{T<:Mixture}(m::T)
  #print("<:Mixture - initializeSampler\n")
  for i in 1:10
    sampleClusters(m)
    sampleAssignments(m)
  end
end

function sample{T<:Mixture}(m::T)
  #print("<:Mixture - sample\n")
  sampleClusters(m)
  sampleAssignments(m)
  sampleNRMIParameters(m)
  sampleClusterHyperparameters(m)
end

function sampleNRMIParameters{T<:Mixture}(m::T)
  #print("<:Mixture - sampleNRMIParameters\n")
  sample(m.nrmi, numData(m), m.clusters)
end

function sampleClusters{T<:Mixture}(m::T)
  #print("<:Mixture - sampleClusters\n")
  for cc in m.clusters
    sample(cc.parameter)
  end
end

function sampleClusterHyperparameters{T<:Mixture}(m::T)
  #print("<:Mixture - sampleClusterHyperparameters\n")
  cdata = Set{Hierarchy}()
  for cc in m.clusters
    push!(cdata, cc.parameter)
  end
  sample(m.prior, cdata)
end

function get{T<:Mixture}(m::T, property::String)
  if property == "nrmi_alpha" return m.nrmi.alpha
  elseif property == "nrmi_sigma" return m.nrmi.sigma
  elseif property == "nrmi_tau" return m.nrmi.tau
  elseif property == "nrmi_logU" return m.nrmi.logU
  elseif property == "prior_precisionInvScale" return m.prior.precisionInvScale
  elseif property == "nbClusters" return length(m.clusters)
  else print("No such property")
  end
end
