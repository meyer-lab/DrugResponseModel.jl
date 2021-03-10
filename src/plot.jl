""" 
        This file contains a function to Plots.plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""
default(size = (900, 400), margin = 0.4cm, legendfontsize = 5, fmt = :pdf)

function G1plots(conc, params, Ylabel; ylim = 2.0)
    conc[1] = 0.05
    concL = log.(conc)
    Plots.plot(concL, params[1, :], xlabel = "log drug conc. [nM]", label = "rate 1", ylabel = Ylabel, lw = 2)
    Plots.plot!(concL, params[2, :], label = "rate 2", lw = 2)
    ylims!(0.0, ylim)
end

function G2plots(conc, params, Ylabel; ylim = 2.0)
    conc[1] = 0.05
    concL = log.(conc)
    Plots.plot(concL, params[1, :], xlabel = "log drug conc. [nM]", label = "rate 1", ylabel = Ylabel, lw = 2)
    Plots.plot!(concL, params[2, :], label = "rate 2", lw = 2)
    Plots.plot!(concL, params[3, :], label = "rate 3", lw = 2)
    Plots.plot!(concL, params[4, :], label = "rate 4", lw = 2)
    ylims!(0.0, ylim)
end


""" Plots.plot the percentage of cells in G2 phase over time for all concentrations on top of each other. Depending on the g2 you pass to it, it could Plots.plot for the data or simulation. """
function plotperc(g2, name, conc, txt)
    time = LinRange(0.0, 95.0, 189)
    p4 = Plots.plot(
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
        Plots.plot!(p4, time, 100.0 .* g2[:, i], lw = 3, color = "black", label = string(conc[i], " nM ", name), alpha = (0.1 * i))
        ylims!((0, 80))
        Plots.plot!(size = (600, 300))
    end
    p4
end

""" Plots.plot the average of three replicates in time series along with the model. it takes in i, that is the concentration number between 1 to 7, j is the number corresponding to the drug number, between 1 to 5. """
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
    Plots.plot!(time, G1[:, i], label = "model G1", color = "darkgreen")
    scatter!(time, g2m[:, i], color = "sienna", label = "data G2", markerstrokewidth = 0.0, markersize = 1.0, alpha = 0.25)
    Plots.plot!(time, G2[:, i], label = "model G2", color = "darkorange")
    ylims!((0.0, 45))
end

function plot_timeCourse(G1, G2, g1m, g2m, i, title)
    x = LinRange(0.0, 95.0, 189)
    pl = Gadfly.Plots.plot(
        layer(x = x, y = G1[:, i], Geom.line, Theme(default_color = colorant"lightblue", line_width = 1.5px)),
        layer(x = x, y = G2[:, i], Geom.line, Theme(default_color = colorant"orange", line_width = 1.5px)),
        layer(x = x, y = g1m[:, i], Geom.line, Theme(default_color = colorant"darkblue", line_width = 1.5px)),
        layer(x = x, y = g2m[:, i], Geom.line, Theme(default_color = colorant"sienna", line_width = 1.5px)),
        layer(x = x, y = G1[:, i] .+ G2[:, i], Geom.line, Theme(default_color = colorant"lightgreen", line_width = 1.5px)),
        layer(x = x, y = g1m[:, i] .+ g2m[:, i], Geom.line, Theme(default_color = colorant"darkgreen", line_width = 1.5px)),
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("time [hr]", orientation = :horizontal),
        Guide.ylabel("cell number", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = 60),
        Guide.manual_color_key("", ["total sim", "total data", "G1 sim", "G1 data", "G2 sim", "G2 data"], ["lightgreen", "darkgreen", "lightblue", "darkblue", "orange", "sienna"]),
        Guide.title(title),
        style(key_position = :inside),
    )
    pl
end
