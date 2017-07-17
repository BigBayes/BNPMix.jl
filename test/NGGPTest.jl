using NRMMM
using Base.Test

alphaShape = 1.0
alphaInvScale = 1.0
sigmaAlpha = 1.0
sigmaBeta = 2.0
tauShape = 1e9
tauInvScale = 1e9

nggp = NGGP(alphaShape, alphaInvScale ,sigmaAlpha, sigmaBeta, tauShape, tauInvScale)

ε = 0.01
@test ≈(logLevy(nggp, 0.5, 5.0), -74.010, atol=ε)
nggp.logU = 5.0
@test logLevy(nggp, 0.5, 0.0) == logLevy(nggp, 0.5)

@test ≈(laplace(nggp, 6.0), 8.225, atol=ε)
nggp.logU = 6.0
@test laplace(nggp, 6.0) == laplace(nggp)

nggp.logU = 4.0
@test ≈(logGamma(nggp, 3, 4.0), -11.116, atol=ε)
@test logGamma(nggp, 3, 4.0) == logGamma(nggp, 3)

@test ≈(logMeanMass(nggp, 7), -2.087, atol=ε)

@test ≈(logMeanTotalMass(nggp, 11), -3.616, atol=ε)

@test ≈(meanNumClusters(nggp, 5), 1.792, atol=ε)
