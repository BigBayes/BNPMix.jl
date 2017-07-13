export
    Cluster,
    isEmpty

mutable struct Cluster#{Theta}
    number::Int
    logmass::Union{Float, Void}
    parameter::Hierarchy#Theta
    w::Float
end

Cluster(number::Int, logmass::Float, parameter::Hierarchy) =
Cluster(number, logmass, parameter, 0.0)

function isEmpty(cc::Cluster)
  return cc.number == 0
end
