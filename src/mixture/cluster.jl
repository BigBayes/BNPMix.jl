export
    Cluster,
    isEmpty

mutable struct Cluster
    number     ::  Int                 # Number of data items assigned to cluster
    logmass    ::  Union{Float, Void}  # Mass (mixture proportion) of cluster
    parameter  ::  Hierarchy           # Parameter
    w          ::  Float               #
    Cluster(number::Int, logmass::Union{Float, Void}, parameter::Hierarchy) =
      new(number, logmass, parameter, 0.0)
end


function isEmpty(cc::Cluster)
  return cc.number == 0
end
