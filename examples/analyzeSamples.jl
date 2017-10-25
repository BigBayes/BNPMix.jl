include("ess.jl")
include("plot.jl")

# Compute statistics
parameters = open(readdlm,"output/nggp.1.6|0.5.galaxy.nc.reuse2")
print("\nESS: ", ess_factor(parameters[:,1])[1])
# parameters2 = open(readdlm,"output/galaxy.nc.slice.parameters")
# print("\nESS: ", ess_factor(parameters2[:,1])[1])

#
plot_histogram(parameters[:,1], "Nb Clusters")
# plot_histograms(parameters[:,6], parameters[:,6],"invscale")
