include("nrmi.jl")

data = [9172.0;9350.0;9483.0;9558.0;9775.0;10227.0;10406.0;16084.0;16170.0;18419.0;18552.0;18600.0;18927.0;19052.0;19070.0;19330.0;19343.0;19349.0;19440.0;19473.0;19529.0;19541.0;19547.0;19663.0;19846.0;19856.0;19863.0;19914.0;19918.0;19973.0;19989.0;20166.0;20175.0;20179.0;20196.0;20215.0;20221.0;20415.0;20629.0;20795.0;20821.0;20846.0;20875.0;20986.0;21137.0;21492.0;21701.0;21814.0;21921.0;21960.0;22185.0;22209.0;22242.0;22249.0;22314.0;22374.0;22495.0;22746.0;22747.0;22888.0;22914.0;23206.0;23241.0;23263.0;23484.0;23538.0;23542.0;23666.0;23706.0;23711.0;24129.0;24285.0;24289.0;24366.0;24717.0;24990.0;25633.0;26960.0;26995.0;32065.0;32789.0;34279.0]
data /= 1e4
data -= mean(data)
mu_0 = 0.0; sigma_1 = 10; sigma_0 = 4*sigma_1;

# alg parameters
alg = "Reuse"
numEmptyClusters = 2

conjugate = false
numBurnin = 5000
numSample = 5000
numThinning = 20

# distributions parameters
alpha = 1.6
sigma = 0.5
tau = 1.0

outputFilename = string("nggp.",alpha,"|",sigma,".galaxy.",conjugate? "cs." : "nc.", lowercase(alg), alg!="Slice"?numEmptyClusters:"")

# initialize prior and factory

prior = NormalGammaIndependent(mu_0, 1/sigma_0^2, 1/sigma_1^2)
factory = NormalNonConjugateFactorySampled()

# build model

nggp = NGGP(alpha, sigma, tau)

if alg == "Neal8" model = MixtureNeal8(nggp, prior, factory, numEmptyClusters)
elseif alg == "Reuse" model = MixtureReuse(nggp, prior, factory, numEmptyClusters)
elseif alg == "Slice" model = MixtureSlice(nggp, prior, factory)
else print("Only Neal8, Reuse and Slice implemented\n")
end

# Add data

numdata = size(data, 1)
newData = zeros(numdata, 1)
newData[:,1] = data
addData(model, newData)

# Parameters STDOUT
print("\nNRMIX ", alg);
print("\n  outputFilename = ", outputFilename)
print("\n  numBurnin = ", numBurnin)
print("\n  numSample = ", numSample)
print("\n  numThinning = ", numThinning)
print("\n  numNewClusters = ", numEmptyClusters)
print("\n  alpha = ", alpha)
print("\n  sigma = ", sigma)
print("\n  tau = ", tau)

# Sampling
collectors = ["nbClusters"]
sampler = Sampler(model, collectors, numBurnin, numSample, numThinning, outputFilename)
runSampler(sampler)
