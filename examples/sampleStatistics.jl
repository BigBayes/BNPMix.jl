using NRMMM

# Compute statistics
parameters = open(readdlm,"output/neal8.parameters")
essNbClusters = ess_factor(parameters[:,1])[1]
print("ESS Julia: ", essNbClusters)

parametersJava = open(readdlm,"../forstefano/galaxy.nc.neal82.parameters")
essNbClustersJava = ess_factor(parametersJava[:,1])[1]
print("\nESS Java: ", essNbClustersJava)

using PlotlyJS
function plot_histograms(x0, x1, title)
    trace1 = histogram(x=x0, opacity=0.75, name="Julia")
    trace2 = histogram(x=x1, opacity=0.75, name="Java")
    data = [trace1, trace2]
    layout = Layout(barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
end
##
plot_histograms(parameters[:,1], parametersJava[:,1], "Nb Clusters")
plot_histograms(parameters[:,2], parametersJava[:,2], "Alpha")
plot_histograms(parametersJava[:,3], parameters[:,3], "Sigma")
#plot_histograms(parametersJava[:,3], parameters[:,3], "Tau")
plot_histograms(parametersJava[:,5], parameters[:,5], "logU")
plot_histograms(parameters[:,6], parametersJava[:,6], "InvScale")
