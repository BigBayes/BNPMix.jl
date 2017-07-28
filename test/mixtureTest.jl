using BNPMix
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

function testMixture(mixtureType)
  model = mixtureType(nggp, prior, factory, numEmptyClusters)

  @test length(model.empties) == numEmptyClusters
  @test length(model.clusters) == 0
  @test model.empties[1].logmass == 0.0

  data = [-1.0;0;1.0;2.0;3.0;4.0;5.0]
  newData = zeros(length(data), 1)
  newData[:,1] = data
  addData(model, newData)

  index = 4
  datum = getDatum(model, index)
  cc = model.map[index]
  @test cc.number == cc.parameter.number
  number = cc.number
  sumX = cc.parameter.sumX
  sumXX = cc.parameter.sumXX
  tt = unassign(model, index, datum)
  @test model.map[index] == nothing
  @test tt == cc
  @test number-1 == tt.number
  @test sumX-datum[1] == tt.parameter.sumX
  @test sumXX-datum[1]*datum[1] == tt.parameter.sumXX

  assign(model, index, datum, cc)
  @test model.map[index] == cc
  @test number == cc.number
  @test sumX == cc.parameter.sumX
  @test sumXX == cc.parameter.sumXX
end

testMixture(MixtureNeal8)
testMixture(MixtureReuse)
