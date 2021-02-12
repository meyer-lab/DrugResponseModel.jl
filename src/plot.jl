""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

default(size = (900, 400), margin = 0.4cm, legendfontsize = 5, fmt = :pdf)

function unit_plot_params(conc, params, stdn, labelY)
    concL = log.(conc)
    plot(
        concL,
        params,
        xlabel = "log drug conc. [nM]",
        label = "",
        ribbon = stdn,
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = labelY,
    )
end

function plot_parameters(concs, parameters, stdn)
    labelYs = [" progression rate [1/hr]", " death rate [1/hr]"]
    pre_labels = ["G1,1", "G1,2", "G2,1", "G2,2", "G2, 3", "G2, 4"]

    p1 = [unit_plot_params(concs, parameters[i, :], stdn[i, :], string(pre_labels[i], labelYs[1])) for i = 1:6]
    p2 = [unit_plot_params(concs, parameters[i, :], stdn[i, :], string(pre_labels[i - 6], labelYs[2])) for i = 7:12]

    plot(p1..., p2..., alpha = 0.6, lw = 2.0, size = (1400, 400), layout = (2, 6), color = [:black :gray])
    ylims!(0.0, 2.0)
end


""" plot the percentage of cells in G2 phase over time for all concentrations on top of each other. Depending on the g2 you pass to it, it could plot for the data or simulation. """
function plotperc(g2, name, conc, txt)
    time = LinRange(0.0, 95.0, 189)
    p4 = plot(
        time,
        100.0 .* g2[:, 1],
        lw = 3,
        color = "red",
        label = "control ",
        legend = :best,
        title = string(txt, " ", name),
        legendfontsize = 5,
        xlabel = "time [hr]",
        ylabel = "S/G2 cell cycle (%) ",
        alpha = 0.7,
    )

    for i = 2:7
        plot!(p4, time, 100.0 .* g2[:, i], lw = 3, color = "black", label = string(conc[i], " nM ", name), alpha = (0.1 * i))
        ylims!((0, 80))
        plot!(size = (600, 300))
    end
    p4
end

""" plot the average of three replicates in time series along with the model. it takes in i, that is the concentration number between 1 to 7, j is the number corresponding to the drug number, between 1 to 5. """
function plotavg(G1, G2, g1m, g2m, i, leg, conc)
    time = LinRange(0.0, 95.0, 189)

    scatter(
        time,
        g1m[:, i],
        title = string(conc, "nM"),
        titlefontsize = 8,
        color = "green",
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = leg,
        label = "data G1",
        xlabel = "time [hr]",
        ylabel = "cell number",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 6,
        alpha = 0.25,
    )
    plot!(time, G1[:, i], label = "model G1", color = "darkgreen")
    scatter!(time, g2m[:, i], color = "sienna", label = "data G2", markerstrokewidth = 0.0, markersize = 1.0, alpha = 0.25)
    plot!(time, G2[:, i], label = "model G2", color = "darkorange")
    ylims!((0.0, 45))
end
