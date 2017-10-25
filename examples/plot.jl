using PlotlyJS

function plot_histogram_with_truth(x0, x1, title="Frequency")
    trace1 = histogram(x=x0, opacity=0.75, name="Empricical")
    trace2 = scatter(;x=1:round(maximum(x0)), y=x1, mode="lines",line_color="rgb(0, 176, 246)",name="MLE fit")
    data = [trace1, trace2]
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
end

function plot_histogram(x0, title="Frequency")
    trace1 = histogram(x=x0, opacity=0.75, name="")
    data = [trace1]
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
end

function plot_histograms(x0, x1, name1="", name2="", title="Frequency")
    trace1 = histogram(x=x0, opacity=0.75, name=name1)
    trace2 = histogram(x=x1, opacity=0.75, name=name2)
    data = [trace1, trace2]
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",barmode="overlay",
    yaxis=attr(title=title,showgrid=false,showline=false),xaxis=attr(title="x",showgrid=false,showline=false))
    plot(data, layout)
end

function linescatter(y0,title="")
    trace1 = scatter(;x=1:length(y0), y=y0, mode="lines+markers")
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",barmode="overlay",yaxis=attr(title=title),xaxis=attr(title="x"))
    plot(data, layout)
    plot([trace1])
end

function linescatter(y0,y1,x0=1:length(y0),x1=1:length(y1),yname="y",xname="x", name1="", name2="")
    trace1 = scatter(;x=x0, y=y0, mode="lines+markers",name=name1)
    trace2 = scatter(;x=x1, y=y1, mode="lines+markers",name=name2)
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",
    yaxis=attr(title=yname,showgrid=false,showline=false),
    xaxis=attr(title=xname,showgrid=false))
    data = [trace1, trace2]
    plot(data, layout)
end

function ESSPlot(y0)
    trace1 = scatter(;x=1:length(y0), y=y0, mode="lines")
    plot([trace1])
    layout = Layout(plot_bgcolor="white",paper_bgcolor="white",
    yaxis=attr(title="Normalized ESS",showgrid=false),xaxis=attr(title="State space time t",showgrid=false))
    plot([trace1], layout)
end

function ESSPlotVariance(Y,name="")
    N = size(Y,2)
    y1 = Y[1,:]
    y2 = Y[2,:]
    y3 = Y[3,:]
    trace2 = scatter(;x=vcat(1:N, N:-1:1),
                     y=vcat(y3, reverse(y1)),
                    #  y=[5.5, 3.0, 5.5, 8.0, 6.0, 3.0, 8.0, 5.0, 6.0, 5.5,
                    #  4.75, 5.0, 4.0, 7.0, 2.0, 4.0, 7.0, 4.4, 2.0, 4.5],
                     fill="tozerox",
                     fillcolor="rgba(0, 176, 246, 0.2)",
                     line_color="transparent",
                     name=name,
                     showlegend=false)

    trace5 = scatter(;x=1:N,
                     y=y2,
                    #  y=[5.0, 2.5, 5.0, 7.5, 5.0, 2.5, 7.5, 4.5, 5.5, 5.],
                     line_color="rgb(0, 176, 246)",
                     mode="lines",
                     name=name)

    data = [trace2, trace5]
    layout = Layout(;paper_bgcolor="white",
                    plot_bgcolor="white",

                    xaxis=attr(title="State space time t",
                               gridcolor="rgb(255, 255, 255)",
                               range=[1, N],
                               showgrid=false,
                               showline=false,
                               showticklabels=true,
                               tickcolor="rgb(127, 127, 127)",
                               ticks="outside",
                               zeroline=false),

                    yaxis=attr(title="Normalized ESS",
                               gridcolor="rgb(255, 255, 255)",
                               showgrid=true,
                               showline=false,
                               showticklabels=true,
                               tickcolor="rgb(127, 127, 127)",
                               ticks="outside",
                               zeroline=false))

    plot(data, layout)
end
