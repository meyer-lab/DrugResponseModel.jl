function plot_pG1_mean_new(y1, y2, y1p, y2p, ymax, Phasename, ylabel, subPlabel, plus)

    # x = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]
    x = [0.25, 0.75, 1.75, 2.25, 3.25, 3.75, 4.75, 5.25, 6.25, 6.75]
    newY1 = zeros(10)
    newY2 = zeros(10)

    for i = 1:5
        newY1[2*i-1] = y1[i] # G1 control
        newY2[2*i-1] = y2[i] # G1 max
        newY1[2*i] = y1p[i] # G2 control
        newY2[2*i] = y2p[i] # G2 max
    end
    conts = scatter(
        x,
        newY1,
        color = "red",
        xlabel = "drugs",
        xrotation = 30,
        label = "Control",
        markerstrokewidth = 0,
        markersize = 8,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$Phasename effects",
    )
    scatter!(
        x,
        newY2,
        color = "cyan4",
        xlabel = "drugs",
        label = "Emax",
        markerstrokewidth = 0,
        markersize = 8,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$Phasename effects",
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.1, ymax))
end
