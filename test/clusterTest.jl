using NRMMM
using Base.Test

prior = NormalGammaIndependent(2.17255, 0.634557, 2.0, 0.2, 6.34557)
factory = NormalNonConjugateFactorySampled()
parameter = construct(factory, prior)
cc = Cluster(0, nothing, parameter, 0.0)

@test isEmpty(cc) == true
cc.number += 1
@test isEmpty(cc) == false
cc.number -= 1
@test isEmpty(cc) == true
