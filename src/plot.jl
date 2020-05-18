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
    ylims!(0.0, 1.2 * maxDeath)

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
    ylims!(0.0, 1.2 * maxDeath)

    plot(p1, p2, p3, p4)
end