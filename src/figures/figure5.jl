""" This plots time series of combination versus their fit. """

p = [0.189371, 1.31225, 0.331636, 3.0, 3.0, 0.274225, 3.32364e-9, 1.48818e-9, 1.5106e-9, 2.04528e-9, 2.3648e-8, 0.0132698, 33.8183, 1.15343, 1.16521e-8, 0.0249639, 0.392687, 0.402105, 2.3868e-8, 0.338073, 0.00245421, 1.36638e-8, 1.62531e-9, 9.78796e-8, 0.0592202, 0.0143732, 0.228744, 2.5, 1.23386, 0.347899, 2.49999, 0.210224]
conc1 = [0, 25, 50, 100, 250];
x1 = ["control" "palbo 50 nM" "palbo 50 nM + lap25 nM" "palbo 50 nM + lap 50 nM" "palbo 50 nM + lap 100 nM" "palbo 50 nM + lap 250 nM"]
x2 = ["control" "palbo 50 nM" "palbo 50 nM + gem 5 nM" "palbo 50 nM + gem 10 nM" "palbo 50 nM + gem 17 nM" "palbo 50 nM + gem 30 nM"]
x3 = ["control" "gem 10 nM" "gem 10 nM + palbo 25 nM" "gem 10 nM + palbo 100 nM" "gem 10 nM + palbo 250 nM" "gem 10 nM + palbo 250 nM"]
x4 = ["control" "gem 10 nM" "gem 10 nM + lap 25 nM" "gem 10 nM + lap 50 nM" "gem 10 nM + lap 100 nM" "gem 10 nM + lap 250 nM"]
x5 = ["control" "dox 20 nM" "dox 20 nM + gem 5 nM" "dox 20 nM + gem 10 nM" "dox 20 nM + gem 17 nM" "dox 20 nM + gem 30 nM"]
x6 = ["control" "lap 100 nM" "lap 100 nM + palbo 25 nM" "lap 100 nM + palbo 100 nM" "lap 100 nM + palbo 250 nM" "lap 100 nM + palbo 250 nM"]
GC = JLD.load("GC.jld")["GC"]
GCsim = zeros(2, 193, 6)

t = LinRange(0.0, 96.0, 193)
combinEffects = DrugResponseModel.my_helper(p, conc1)
for i = 1:6 # concentration number
    GCsim[1, :, i], GCsim[2, :, i], _ = predict(combinEffects[:, i], combinEffects[:, 1], t)
end

function plot_fig5(g1, g1data, labels, tite, G, subPlabel)
    time = LinRange(0.0, 96.0, 193)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = labels,
        fg_legend = :transparent,
        palette = :PuBu_6,
        title = tite,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        xlabel = "time [hr]",
        xticks = 0:24.0:96.0,
        ylabel = "$G cell number",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(time, g1data, lw = 2, linestyle = :dot, label = ["" "" "" "" "" "" ""])
    annotate!(-25.0, 57.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 50))
    p
end

function plot_pCombo(efcs, ymax, label2, Phasename, ylabel, subPlabel, plus)

    x = ["G11", "G12", "G21", "G22", "G23", "G24"]
    y1 = efcs[1:6, 1] # control
    y2 = efcs[1:6, 2] # base drug alone
    y3 = efcs[1:6, 6] # combination at max of drug B
    scatter(
        x,
        y1,
        color = "blue",
        xlabel = "sub-phase",
        xrotation = 40,
        label = "Control",
        markerstrokewidth = 0,
        markershape=:circle,
        markersize = 10,
        alpha=0.6,
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
        y2,
        color = "red",
        alpha=0.6,
        label = label2,
        markerstrokewidth = 0,
        markershape=:square,
        markersize = 8,
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
    )
    scatter!(
        x,
        y3,
        color = "green",
        label = "Emax Combo",
        alpha=0.6,
        markerstrokewidth = 0,
        markershape=:star5,
        markersize = 8,
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.1, ymax))
end

function figure5()
    p1 = plot_fig5(GCsim[1, :, :], GC[1, :, :, 1], x1, "Palbo 50 - Lap combo", "G1", "A")
    p2 = plot_fig5(GCsim[2, :, :], GC[2, :, :, 1], x1, "Palbo 50 - Lap combo", "G2", "B")
    p3 = plot_pCombo(combinEffects[1:6, :], 3.75, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "C", 0.35)
    p4 = plot_pCombo(combinEffects[7:12, :], 0.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "D", 0.06)

    fig5 = plot(p1, p2, p3, p4, size = (1050, 800), layout = (2, 2))
    savefig(fig5, "figure5.svg")
end