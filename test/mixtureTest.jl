using NRMMM
using Base.Test

prior = NormalGammaIndependent(2.17255, 0.634557, 2.0, 0.2, 6.34557)
factory = NormalNonConjugateFactorySampled()
parameter = construct(factory, prior)

alphaShape = 1.0
alphaInvScale = 1.0
sigmaAlpha = 1.0
sigmaBeta = 2.0
tauShape = 1e9
tauInvScale = 1e9

nggp = NGGP(alphaShape, alphaInvScale ,sigmaAlpha, sigmaBeta, tauShape, tauInvScale)

numEmptyClusters = 2
model = MixtureNeal8(nggp, prior, factory, numEmptyClusters)

@test length(model.empties) == numEmptyClusters
@test length(model.clusters) == 0
@test model.empties[1].logmass == 0.0

#data = [-1.0;0;1.0;2.0]
