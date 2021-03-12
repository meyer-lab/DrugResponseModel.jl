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

function plot_G1_params(pr, i, titles, ymax)

    # progressions = [12 x 8]
    G1p = vcat(pr[i:i+1, 1, 1], pr[i:i+1, 1, 2], pr[i:i+1, 1, 3], pr[i:i+1, 1, 4], pr[i:i+1, 1, 5]); # G1 progression control
    G1p2 = vcat(pr[i:i+1, 8, 1], pr[i:i+1, 8, 2], pr[i:i+1, 8, 3], pr[i:i+1, 8, 4], pr[i:i+1, 8, 5]); # G1 progression max

    g1 = repeat(["G11", "G12"], 5)
    xx = zeros(10)
    for i=1:5
        xx[(2*i-1)] = i - 0.15
        xx[2*i] = i + 0.15
    end
    xticks = [1.0, 2.0, 3.0, 4.0, 5.0]
    labels = Dict(zip(xticks, ["Lpt","Dox","Gemc","Taxol", "palbo"]))
    df1 = DataFrame(x=xx, ymin=G1p, x2=xx.+0.05, ymax=G1p2, subphase=g1)
    # df = vcat(df1, df2)
    pl = Gadfly.plot(df1, xmin=:x, ymin=:ymin, xmax=:x2, ymax=:ymax, color=:subphase, Geom.rect,
        Stat.xticks(ticks = xticks),
        Scale.x_continuous(labels = x -> labels[x]),
        Guide.xlabel("Drugs", orientation = :horizontal),
        Guide.ylabel("rates 1/[hr]", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = ymax),
        Guide.title(titles),
        style(key_position = :top),
    )
    pl
end

function plot_G2_params(pr, i, titles, ymax)

    # progressions = [12 x 8]
    G2p = vcat(pr[i:i+3, 1, 1], pr[i:i+3, 1, 2], pr[i:i+3, 1, 3], pr[i:i+3, 1, 4], pr[i:i+3, 1, 5]); # G2 death control
    G2p2 = vcat(pr[i:i+3, 8, 1], pr[i:i+3, 8, 2], pr[i:i+3, 8, 3], pr[i:i+3, 8, 4], pr[i:i+3, 8, 5]); # G2 death max
    g2 = repeat(["G21", "G22", "G23", "G24"], 5)
    xx = [0.25, 0.75, 1.25, 1.75]
    xs = vcat(xx, xx .+ 4, xx .+ 8, xx .+ 12, xx .+ 16)

    xticks = [1.0, 5.0, 9.0, 13.0, 17.0]
    labels = Dict(zip(xticks, ["Lpt","Dox","Gemc","Taxol", "palbo"]))
    df1 = DataFrame(x=xs, ymin=G2p, x2=xs.+0.2, ymax=G2p2, subphase=g2)
    pl = Gadfly.plot(df1, xmin=:x, ymin=:ymin, xmax=:x2, ymax=:ymax, color=:subphase, Geom.rect,
        Stat.xticks(ticks = xticks),
        Scale.x_continuous(labels = x -> labels[x]),
        Guide.xlabel("Drugs", orientation = :horizontal),
        Guide.ylabel("rates 1/[hr]", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = ymax),
        Guide.title(titles),
        style(key_position = :top),
    )
    pl
end