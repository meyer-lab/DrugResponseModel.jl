""" This plots time series of combination versus their fit. """

p = [0.183721, 0.193429, 2.10987, 0.492517, 2.26232, 0.0630642, 0.000118551, 2.71169e-5, 0.000224803, 7.33529e-5, 0.000683488, 1.34972e-5, 0.467268, 0.285912, 2.15848, 0.466075, 0.316368, 0.580013, 0.000562746, 4.35196e-5, 0.000121694, 3.02369e-5, 1.49331e-5, 0.000504069, 78.9947, 1.07977, 2.48639, 0.00865033, 1.25326, 2.48478, 0.152791, 1.39507, 0.00151377, 0.00225855, 0.000307931, 0.129873, 0.0120839, 0.000253296, 8.50684, 5.43357, 1.96648, 0.411635, 0.338033, 1.4666, 0.285535, 0.899417, 0.150104, 0.00175209, 0.00275891, 0.013657, 0.0125972, 0.000519278, 5.50678, 0.266267, 2.48314, 0.0065199, 0.167905, 0.0480877, 0.255671, 2.45708, 0.0156757, 4.98363e-5, 1.20754e-5, 3.99582e-5, 1.5227e-5, 0.0640394, 0.519009, 0.368302, 0.466753, 2.29984, 0.248386, 0.440042];

concs = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
x1 = ["control" "palbo 50 nM" "palbo 50 nM + lap25 nM" "palbo 50 nM + lap 50 nM" "palbo 50 nM + lap 100 nM" "palbo 50 nM + lap 250 nM"]
x2 = ["control" "palbo 50 nM" "palbo 50 nM + gem 5 nM" "palbo 50 nM + gem 10 nM" "palbo 50 nM + gem 17 nM" "palbo 50 nM + gem 30 nM"]
x3 = ["control" "gem 10 nM" "gem 10 nM + palbo 25 nM" "gem 10 nM + palbo 50 nM" "gem 10 nM + palbo 100 nM" "gem 10 nM + palbo 250 nM"]
x4 = ["control" "gem 10 nM" "gem 10 nM + lap 25 nM" "gem 10 nM + lap 50 nM" "gem 10 nM + lap 100 nM" "gem 10 nM + lap 250 nM"]
x5 = ["control" "lap 100 nM" "lap 100 nM + palbo 25 nM" "lap 100 nM + palbo 50 nM" "lap 100 nM + palbo 100 nM" "lap 100 nM + palbo 250 nM"]
x6 = ["control" "lap 100 nM" "lap 100 nM + gem 5 nM" "lap 100 nM + gem 10 nM" "lap 100 nM + gem 17 nM" "lap 100 nM + gem 30 nM"]
x7 = ["control" "dox 20 nM" "dox 20 nM + gem 5 nM" "dox 20 nM + gem 10 nM" "dox 20 nM + gem 17 nM" "dox 20 nM + gem 30 nM"]
x8 = ["control" "Pax 2 nM" "pax 2 nM + palbo 50 nM" "pax 2 nM + dox 20 nM" "pax 2 nM + lap 50" "pax 2 nM + lap100" "pax 2 nM + gem10 nM"]

g = JLD.load("g.jld")["g"]
gT = g[1, :, :, :] .+ g[2, :, :, :]
GCsim = zeros(2, 193, 6, 8)

t = LinRange(0.0, 96.0, 193)

combinEffects = DrugResponseModel.my_helper(p)
single_effs = zeros(12, 5, 5)
single_effs[:, 2, 2] = p[1:12]
single_effs[:, 2, 4] = p[13:24]
single_effs[:, :, [1, 3, 5]] = getODEparams(p[25:end], concs);

for j=1:8
    for i = 1:6 # concentration number
        GCsim[1, :, i, j], GCsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
    end
end
GCsimT = GCsim[1, :, :, :] .+ GCsim[2, :, :, :] # Total

function plot_fig5(time, g1, g1data, labels, tite, G, subPlabel, palet)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = labels,
        fg_legend = :transparent,
        palette = palet,
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
# p = [0.621894, 0.327901, 1.43602, 0.138798, 2.79013, 1.32374, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

function plot_pCombo(efcs, single_effs, drugB, ymax, label2, Phasename, ylabel, subPlabel, plus)

    x = ["G11", "G12", "G21", "G22", "G23", "G24"]
    y1 = efcs[1:6, 1] # control
    y2 = efcs[1:6, 2] # base drug alone
    y3 = efcs[1:6, 5] # combination at max of drug B 
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
        alpha=0.4,
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
        alpha=0.4,
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
        alpha=0.4,
        markerstrokewidth = 0,
        markershape=:star5,
        markersize = 8,
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
    )
    scatter!(
        x,
        single_effs,
        color = "black",
        label = "Emax $drugB",
        alpha=0.6,
        markerstrokewidth = 0,
        markershape=:+,
        markersize = 8,
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.1, ymax))
end

""" Plot G1 and S/G2 cell numbers of experiments vs simulations. """
function figure5()
    p1 = plot_fig5(t, GCsim[1, :, :, 1], g[1, :, :, 1], x1, "Palbo 50 + Laps", "G1", "A", :PuBu_6)
    p2 = plot_fig5(t, GCsim[2, :, :, 1], g[2, :, :, 1], x1, "Palbo 50 + Laps", "G2", "B", :PuBu_6)
    p3 = plot_fig5(t, GCsim[1, :, :, 2], g[1, :, :, 2], x2, "Palbo 50 + Gems", "G1", "C", :PuBu_6)
    p4 = plot_fig5(t, GCsim[2, :, :, 2], g[2, :, :, 2], x2, "Palbo 50 + Gems", "G2", "D", :PuBu_6)
    p5 = plot_fig5(t, GCsim[1, :, :, 3], g[1, :, :, 3], x3, "Gem 10 + Palbos", "G1", "E", :PuBu_6)
    p6 = plot_fig5(t, GCsim[2, :, :, 3], g[2, :, :, 3], x3, "Gem 10 + Palbos", "G2", "F", :PuBu_6)
    p7 = plot_fig5(t, GCsim[1, :, :, 4], g[1, :, :, 4], x4, "Gem 10 + Laps", "G1", "G", :PuBu_6)
    p8 = plot_fig5(t, GCsim[2, :, :, 4], g[2, :, :, 4], x4, "Gem 10 + Laps", "G2", "H", :PuBu_6)
    p9 = plot_fig5(t, GCsim[1, :, :, 5], g[1, :, :, 5], x5, "Lap 100 + Palbos", "G1", "I", :PuBu_6)
    p10 = plot_fig5(t, GCsim[2, :, :, 5], g[2, :, :, 5], x5, "Lap 100 + Palbos", "G2", "J", :PuBu_6)
    p11 = plot_fig5(t, GCsim[1, :, :, 6], g[1, :, :, 6], x6, "Lap 100 + Gems", "G1", "K", :PuBu_6)
    p12 = plot_fig5(t, GCsim[2, :, :, 6], g[2, :, :, 6], x6, "Lap 100 + Gems", "G2", "L", :PuBu_6)
    p13 = plot_fig5(t, GCsim[1, :, :, 7], g[1, :, :, 7], x7, "Dox 20 + Gem", "G1", "M", :PuBu_6)
    p14 = plot_fig5(t, GCsim[2, :, :, 7], g[2, :, :, 7], x7, "Dox 20 + Gem", "G2", "N", :PuBu_6)
    p15 = plot_fig5(t, GCsim[1, :, :, 8], g[1, :, :, 8], x8, "Pax 2 + others", "G1", "o", :PuBu_6)
    p16 = plot_fig5(t, GCsim[2, :, :, 8], g[2, :, :, 8], x8, "Pax 2 + others", "G2", "P", :PuBu_6)

    fig5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, size = (1900, 1600), layout = (4, 4))
    savefig(fig5, "figure5.svg")
end

""" Plot Total cell numbers of experiments vs simulation. """
function figure5T()
    p1 = plot_fig5(t, GCsimT[:, :, 1], gT[:, :, 1], x1, "Palbo 50 + Laps", "Total", "A", :PuBu_6)
    p2 = plot_fig5(t, GCsimT[:, :, 2], gT[:, :, 2], x2, "Palbo 50 + Gems", "Total", "B", :PuBu_6)
    p3 = plot_fig5(t, GCsimT[:, :, 3], gT[:, :, 3], x3, "Gem 10 + Palbos", "Total", "C", :PuBu_6)
    p4 = plot_fig5(t, GCsimT[:, :, 4], gT[:, :, 4], x4, "Gem 10 + Laps", "Total", "D", :PuBu_6)
    p5 = plot_fig5(t, GCsimT[:, :, 5], gT[:, :, 5], x5, "Lap 100 + Palbos", "Total", "E", :PuBu_6)
    p6 = plot_fig5(t, GCsimT[:, :, 6], gT[:, :, 6], x6, "Lap 100 + Gems", "Total", "F", :PuBu_6)
    p7 = plot_fig5(t, GCsimT[:, :, 7], gT[:, :, 7], x7, "Dox 20 + Gem", "Total", "G", :PuBu_6)
    p8 = plot_fig5(t, GCsimT[:, :, 8], gT[:, :, 8], x8, "Pax 2 + others", "Total", "H", :PuBu_6)

    fig5t = plot(p1, p2, p3, p4, p5, p6, p7, p8, size = (1800, 1000), layout = (2, 4))
    savefig(fig5t, "figure5T.svg")
end

""" Plot inferred rate parameters. """
function figureS5()
    p1 = plot_pCombo(combinEffects[1:6, :, 1], single_effs[1:6, end, 1], "Lap", 10.95, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "A", 0.35)
    p2 = plot_pCombo(combinEffects[7:12, :, 1], single_effs[1:6, end, 1], "Lap", 2.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "B", 0.06)
    p3 = plot_pCombo(combinEffects[1:6, :, 2], single_effs[1:6, end, 3], "Gem", 10.95, "Palbo 50 nM", "Palbo 50 nM + Gem", "progression rates [1/hr]", "C", 0.35)
    p4 = plot_pCombo(combinEffects[7:12, :, 2], single_effs[1:6, end, 3], "Gem", 2.5, "Palbo 50 nM", "Palbo 50 nM + Gem", "death rates [1/hr]", "D", 0.06)
    p5 = plot_pCombo(combinEffects[1:6, :, 3], single_effs[1:6, end, 5], "Palbo", 10.95, "Gem 10 nM", "Gem 10 nM + Palbo", "progression rates [1/hr]", "E", 0.35)
    p6 = plot_pCombo(combinEffects[7:12, :, 3], single_effs[1:6, end, 5], "Palbo", 2.5, "Gem 10 nM", "Gem 10 nM + Palbo", "death rates [1/hr]", "F", 0.06)
    p7 = plot_pCombo(combinEffects[1:6, :, 4], single_effs[1:6, end, 1], "Lap", 10.95, "Gem 10 nM", "Gem 10 nM + Lap", "progression rates [1/hr]", "G", 0.35)
    p8 = plot_pCombo(combinEffects[7:12, :, 4], single_effs[1:6, end, 1], "Lap", 2.5, "Gem 10 nM", "Gem 10 nM + Lap", "death rates [1/hr]", "H", 0.06)
    p9 = plot_pCombo(combinEffects[1:6, :, 5], single_effs[1:6, end, 1], "Lap", 10.95, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "I", 0.35)
    p10 = plot_pCombo(combinEffects[7:12, :, 5], single_effs[7:12, end, 1], "Lap", 2.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "J", 0.06)
    p11 = plot_pCombo(combinEffects[1:6, :, 6], single_effs[1:6, end, 3], "Gem", 10.95, "Lap 100 nM", "Lap 100 nM + Gem", "progression rates [1/hr]", "K", 0.35)
    p12 = plot_pCombo(combinEffects[7:12, :, 6], single_effs[7:12, end, 3], "Gem", 2.5, "Lap 100 nM", "Lap 100 nM + Gem", "death rates [1/hr]", "L", 0.06)
    p13 = plot_pCombo(combinEffects[1:6, :, 7], single_effs[1:6, end, 3], "Gem", 10.95, "Dox 20 nM", "Dox 20 nM + Gem", "progression rates [1/hr]", "M", 0.35)
    p14 = plot_pCombo(combinEffects[7:12, :, 7], single_effs[7:12, end, 3], "Gem", 2.5, "Dox 20 nM", "Dox 20 nM + Gem", "death rates [1/hr]", "N", 0.06)
    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)

    figS5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p0, p0, size = (1800, 1600), layout = (4, 4))
    savefig(figS5, "figureS5.svg")
end