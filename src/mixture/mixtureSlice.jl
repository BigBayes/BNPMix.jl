mutable struct MixtureSlice <: Mixture
    nrmi::NRMI
    prior::Prior
    factory::Factory
    numNewClusters::Int
    data::Union{Void, Array{Float, 2}}
    map::Union{Array{Union{Cluster, Void}, 0}, Array{Union{Cluster, Void}, 1}}
    clusters::Set{Cluster}
    numClustersRatio::Float
    minslice::Float
    slice::Union{Array{Float, 0},Array{Float, 1}}
    clusterarray::Union{Array{Float, 0},Array{Float, 1}}
    # Constructor
    function MixtureReuse(nrmi::NRMI, prior::Prior, factory::Factory, numNewClusters::Int)
      this = new(nrmi, prior, factory, numNewClusters, nothing, Array{Union{Void, Cluster}}(0), Set{Cluster}(), 5.0, 1.0e-8, Array{Float}(0), Array{Float}(0))
      return this
    end
end

function addData(m::MixtureReuse, traindata::Array{Float, 2})
  #addDataAbstract(m, traindata)
  #sampleAssignments(m)
end

function sampleAssignments(m::MixtureReuse)
  numData = size(m.data, 1)
  for i in 1:numData
    sampleAssignment(m, i)
  end
  m.clusters = Set()
  num = 0
  for cc in m.clusterarray
    if !isEmpty(cc)
      push!(clusters, cc)
    end #destruct(factory, prior, cc.parameter)
    num += cc.number
  end
  assert(num == numData)
end

function sampleAssignment(m::MixtureReuse, index::Int)
  datum = getDatum(m, index)
  unassign(m, index, datum)
  thisslice = slice[index]

  for cc in m.clusterarray
    if cc.logmass < thisslice break end
    cc.w = logPredictive(cc.parameter, datum)
  end

  ma = -Inf
  for cc in m.clusterarray
    if cc.logmass < thisslice break end
    if ma < cc.w ma = cc.w
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
