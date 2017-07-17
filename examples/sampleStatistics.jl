using NRMMM

# Compute statistics
parameters = open(readdlm,"/Users/EmileMathieu/code/forstefano/galaxy.nc.reuse.parameters")
print("\nESS: ", ess_factor(parameters[:,1])[1])

parameters2 = open(readdlm,"/Users/EmileMathieu/code/forstefano/galaxy.nc.neal8.parameters")
print("\nESS: ", ess_factor(parameters2[:,1])[1])

parameters3 = open(readdlm,"/Users/EmileMathieu/code/forstefano/galaxy.nc.slice.parameters")
print("\nESS: ", ess_factor(parameters3[:,1])[1])

parameters4 = open(readdlm,"output/galaxy.nc.slice.parameters")
parameters5 = open(readdlm,"output/galaxy.nc.reuse2.parameters")

using PlotlyJS
function plot_histograms(x0, x1, title)
    trace1 = histogram(x=x0, opacity=0.75, name="1")
    trace2 = histogram(x=x1, opacity=0.75, name="2")
    data = [trace1, trace2]
    layout = Layout(barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
end
##
plot_histograms(parameters4[:,1], parameters5[:,1], "Nb Clusters")
plot_histograms(parameters4[:,2], parameters5[:,2], "Alpha")
plot_histograms(parameters4[:,3], parameters5[:,3], "Sigma")
#plot_histograms(parameters4[:,3], parameters5[:,3], "Tau")
plot_histograms(parameters4[:,5], parameters5[:,5], "logU")
plot_histograms(parameters4[:,6], parameters5[:,6], "InvScale")
