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

function plot_percTimeCourse(G1, G2, g1m_min, g1m_max, concs, title)
    x = LinRange(0.0, 95.0, 189)

    G1_perc = (100 .* G1) ./ (G1 .+ G2)
    df1 = DataFrame(x=x, y=G1_perc[:, 1], ymin=g1m_min[:, 1], ymax=g1m_max[:, 1], f="control")
    df2 = DataFrame(x=x, y=G1_perc[:, 3], ymin=g1m_min[:, 3], ymax=g1m_max[:, 3], f="$(concs[3]) nM")
    df3 = DataFrame(x=x, y=G1_perc[:, 5], ymin=g1m_min[:, 5], ymax=g1m_max[:, 5], f="$(concs[5]) nM")
    df4 = DataFrame(x=x, y=G1_perc[:, 7], ymin=g1m_min[:, 7], ymax=g1m_max[:, 7], f="$(concs[7]) nM")

    df = vcat(df1, df2, df3, df4)
    pl = Gadfly.plot(df, x=:x, y=:y, ymin=:ymin, ymax=:ymax, color=:f, Geom.line, Geom.ribbon,
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("time [hr]", orientation = :horizontal),
        Guide.ylabel("G1 %", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = 100),
        Guide.title(title),
        style(key_position = :inside),
    )
    pl
end

function plot_totalTimeCourse(total_sim, tots_min, tots_max, concs, title)
    x = LinRange(0.0, 95.0, 189)

    df1 = DataFrame(x=x, y=total_sim[:, 1], ymin=tots_min[:, 1], ymax=tots_max[:, 1], f="control")
    df2 = DataFrame(x=x, y=total_sim[:, 3], ymin=tots_min[:, 3], ymax=tots_max[:, 3], f="$(concs[3]) nM")
    df3 = DataFrame(x=x, y=total_sim[:, 5], ymin=tots_min[:, 5], ymax=tots_max[:, 5], f="$(concs[5]) nM")
    df4 = DataFrame(x=x, y=total_sim[:, 7], ymin=tots_min[:, 7], ymax=tots_max[:, 7], f="$(concs[7]) nM")

    df = vcat(df1, df2, df3, df4)
    pl = Gadfly.plot(df, x=:x, y=:y, ymin=:ymin, ymax=:ymax, color=:f, Geom.line, Geom.ribbon,
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("time [hr]", orientation = :horizontal),
        Guide.ylabel("cell number", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = 60),
        Guide.title(title),
        style(key_position = :inside),
    )
    pl
end

function plot_progression_rates(progressions, concs, title)

    df1 = DataFrame(x=concs, y=progressions, f="p1")
    df2 = DataFrame(x=concs, y=progressions, f="p2")
    df3 = DataFrame(x=concs, y=progressions, f="p3")
    df4 = DataFrame(x=concs, y=progressions, f="p4")

    df = vcat(df1, df2, df3, df4)
    pl = Gadfly.plot(df, x=:x, y=:y, color=:f, Geom.line,
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("concentrations [hr]", orientation = :horizontal),
        Guide.ylabel("cell number", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = 3),
        Guide.title(title),
        style(key_position = :inside),
    )
    pl
end