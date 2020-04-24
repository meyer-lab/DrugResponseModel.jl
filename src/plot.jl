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
    ylims!(0.0, 1.2 * maximum(parameters[1, :]))

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
    ylims!(0.0, 1.2 * maximum(parameters[2, :]))

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

function plot2(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3, g1s1, g1s2, g1s3, g2s1, g2s2, g2s3, i, j)
    time = LinRange(0.0, 95.0, 189)
    meang1 = zeros(189, 8, 5)
    meang2 = zeros(189, 8, 5)
    stdg1 = zeros(189, 8, 5)
    stdg2 = zeros(189, 8, 5)
    for k = 1:5
        meang1[:, :, k], meang2[:, :, k], stdg1[:, :, k], stdg2[:, :, k] =
            mean_std_data(g1s1[:, :, k], g1s2[:, :, k], g1s3[:, :, k], g2s1[:, :, k], g2s2[:, :, k], g2s3[:, :, k])
    end
    plot(time, meang1[:, i, j]; ribbon = stdg1[:, i, j], color = 1, label = "", xlabel = "time [hr]", ylabel = "cell number", alpha = 0.1)
    plot!(time, G1_1, label = "G1", color = 1)
    plot!(time, G1_2, label = "", color = 1)
    plot!(time, G1_3, label = "", color = 1)
    plot!(time, meang2[:, i, j]; ribbon = stdg2[:, i, j], color = 2, label = "", alpha = 0.1)
    plot!(time, G2_1, label = "G2", color = 2)
    plot!(time, G2_2, label = "", color = 2)
    plot!(time, G2_3, label = "", color = 2)
    ylims!((0.0, 45))
end
