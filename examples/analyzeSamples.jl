using BNPMix
using PlotlyJS
include("ess.jl")

# Compute statistics
parameters = open(readdlm,"output/galaxy.nc.reuse2.parameters")
print("\nESS: ", ess_factor(parameters[:,1])[1])
parameters2 = open(readdlm,"output/galaxy.nc.slice.parameters")
print("\nESS: ", ess_factor(parameters2[:,1])[1])

function plot_histograms(x0, x1, title)
    trace1 = histogram(x=x0, opacity=0.75, name="1")
    trace2 = histogram(x=x1, opacity=0.75, name="2")
    data = [trace1, trace2]
    layout = Layout(barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
end

plot_histograms(parameters[:,1], parameters2[:,1], "Nb Clusters")
