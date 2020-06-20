""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

default(size = (900, 400), margin = 0.4cm, legendfontsize = 7, fmt = :pdf)

#-------------------------- plot for the G1 G2 correlation ---------------------------#
function correlationPlot(g1::Matrix, g2::Matrix, labels::Array, xlabel::String, ylabel::String, ymax::Int)
    pl = [
        scatter(
            g1[:, i],
            g2[:, i],
            title = labels[i],
            xlabel = xlabel,
            ylabel = ylabel,
            alpha = 0.6,
            color = [:gray],
            line = :dot,
            marker = (:dot, 1.5),
        ) for i = 1:8
    ]
    plot(pl..., legend = :false, layout = (2, 4))
    ylims!((0, ymax))
    xlims!((0, ymax))
end

function plot_parameters(conc_l, parameters, stdn)
    conc = log.(conc_l)
    maxProg = maximum(parameters[1:2, :])
    p1 = plot(
        conc,
        parameters[1, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        ribbon = stdn[1, :],
        lw = 2.0,
        alpha = 0.6,
        color = [:black :gray],
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "G1 progression rate",
    )
    ylims!(0.0, 1.2 * maxProg)

    p2 = plot(
        conc,
        parameters[2, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        ribbon = stdn[2, :],
        lw = 2.0,
        alpha = 0.6,
        color = [:black :gray],
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "G2 progression rate",
    )
    ylims!(0.0, 1.2 * maxProg)

    maxDeath = maximum(parameters[3:4, :])

    p3 = plot(
        conc,
        parameters[3, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        ribbon = stdn[3, :],
        lw = 2.0,
        alpha = 0.6,
        color = [:black :gray],
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "G1 death rate",
    )
    ylims!(0.0, 0.8)

    p4 = plot(
        conc,
        parameters[4, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        ribbon = stdn[4, :],
        lw = 2.0,
        alpha = 0.6,
        color = [:black :gray],
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "G2 death rate",
    )
    ylims!(0.0, 0.8)
    plot(p1, p2, p3, p4)
end

""" plot the percentage of cells in G2 phase over time for all concentrations on top of each other. Depending on the g2 you pass to it, it could plot for the data or simulation. """
function plotperc(g2, name, conc)
    time = LinRange(0.0, 95.0, 189)
    p4 = plot(time, 100.0 .* g2[:, 1], lw=3, color="red", label = "control ", legend=:best, title=string(name), legendfontsize=5, xlabel = "time [hr]", ylabel = "S/G2 cell cycle (%) ", alpha = 0.7)

    for i=2:7
        plot!(p4, time, 100.0 .* g2[:, i], lw=3, color="black", label = string(conc[i], " nM ", name), alpha = (1.0 - 0.1*i))
        ylims!((0,80))
        plot!(size=(300, 300), dpi=150)
    end
    p4
end

""" plot the average of three replicates in time series along with the model. it takes in i, that is the concentration number between 1 to 7, j is the number corresponding to the drug number, between 1 to 5. """
function plotavg(G1, G2, g1m, g2m, i, leg, conc)
    time = LinRange(0.0, 95.0, 189)

    scatter(time, g1m[:, i], title=string(conc, "nM"), color = "green", markersize=1.0, markerstrokewidth=0, legend=leg, label = "data G1", xlabel = "time [hr]", ylabel = "cell number", alpha = 0.8)
    plot!(time, G1[:, i], label = "model G1", color = "darkgreen")
    scatter!(time, g2m[:, i], color = "sienna", label = "data G1", markerstrokewidth=0, markersize=1.0, alpha = 0.8)
    plot!(time, G2[:, i], label = "model G2", color = "darkorange")
    ylims!((0.0, 45))
end
