include("../xfamily/normalGammaIndependent.jl")

mutable struct Cluster#{Theta}
    number::Int64
    logmass::Union{Float64, Void}
    parameter::Hierarchy#Theta
    w::Float64
end

Cluster(number::Int64, logmass::Float64, parameter::Hierarchy) =
Cluster(number, logmass, parameter, 0.0)

function isEmpty(cc::Cluster)
  return cc.number == 0
end
