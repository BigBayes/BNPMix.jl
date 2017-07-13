include("mixture.jl")
include("../nrmi/NGGP.jl")

mutable struct MixtureNeal8 <: Mixture
    nrmi::NRMI
    prior::Prior
    factory::Factory
    numNewClusters::Int64
    data::Union{Array{Float64, 0}, Array{Float64, 1}, Array{Float64, 2}}
    map::Union{Array{Union{Cluster, Void}, 0}, Array{Union{Cluster, Void}, 1}}
    clusters::Set{Cluster}
    numClustersRatio::Float64
    empties::Array{Union{Cluster, Void}, 1}

    # Constructor
    function MixtureNeal8(nrmi::NRMI, prior::Prior, factory::Factory, numNewClusters::Int64)
      this = new(nrmi, prior, factory, numNewClusters, zeros(Float64, 0), Array(Union{Void, Cluster}, 0), Set(), 5.0, Array(Union{Void, Cluster}, 0))
      initializeEmpties(this)
      return this
    end
end

function initializeEmpties(m::MixtureNeal8)
  for i in 1:m.numNewClusters
    cc = newCluster(m, 0.0)
    push!(m.empties, cc)
  end
end

function addData(m::MixtureNeal8, traindata::Union{Array{Float64, 1}, Array{Float64, 2}})
  addDataAbstract(m, traindata)
  #sampleAssignments(m)
end

function sampleAssignments(m::MixtureNeal8)
  #print("MixtureNeal8 - Sample Assignments\n")
  numData = size(m.data, 1)
  for cc in m.clusters
    cc.logmass = logMeanMass(m.nrmi, cc.number)
  end
  for i in 1:m.numNewClusters
    cc = m.empties[i]
    cc.logmass = nothing
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

function sampleAssignment(m::MixtureNeal8, index::Int64)
  datum = getDatum(m, index)
  tt = unassign(m, index, datum)

  # initialize empty clusters
  new0 = m.empties[1]
  if !isEmpty(tt)
    tt.logmass = logMeanMass(m.nrmi, tt.number)
    sample(new0.parameter)
  else
    delete!(m.clusters, tt)
    #factory.destruct(prior, new0.parameter) # useless
    new0.parameter = tt.parameter
  end

  for i in 1:m.numNewClusters
    sample(m.empties[i].parameter)
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
