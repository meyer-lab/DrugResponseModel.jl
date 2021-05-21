""" This plots time series of combination versus their fit. """

# p= [0.143126, 1.56756, 3.16138, 0.378522, 1.66322, 0.0920677, 0.000339941, 0.000466815, 0.000498842, 0.000378931, 0.00024666, 2.87482e-5, 0.26904, 2.21656, 1.00714, 3.37664, 1.85998, 0.209481, 1.71511e-5, 0.000301283, 0.000109804, 0.112384, 0.00076013, 0.00199649, 41.3821, 1.73555, 0.0633362, 0.0623595, 1.13678, 0.000480683, 1.05972, 0.358358, 0.000934846, 0.00453242, 3.43235e-5, 0.000236285, 0.0104269, 0.0168968, 6.97406, 5.27271, 0.803043, 0.677746, 0.783597, 4.65155, 0.637839, 0.321503, 0.000146549, 8.4454e-5, 0.000555598, 0.0169459, 2.77595e-5, 0.0370271, 5.05281, 0.669549, 0.168418, 0.264814, 0.206935, 4.74246, 0.533856, 0.302541, 2.20575e-5, 6.81311e-5, 0.000175995, 0.0603435, 0.000106947, 0.0109141, 0.226633, 3.91577, 0.906695, 0.337063, 4.99587, 0.214753]
p = [0.547678, 1.56455, 0.470488, 1.04509, 2.0, 0.449469, 0.01, 0.0100001, 0.01, 0.0272842, 0.0100002, 0.01, 0.46596, 1.64378, 1.42913, 0.809924, 1.48764, 0.534908, 0.01, 0.01, 0.0389332, 0.01, 0.01, 0.01, 99.8637, 1.3282, 0.01, 0.0100001, 0.0215467, 1.35121, 1.8, 0.358891, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 16.2989, 1.2226, 0.2905, 0.01, 0.194664, 0.0100005, 0.934234, 0.410242, 0.01, 0.01, 0.01, 0.01, 0.0100001, 0.01, 390.204, 0.400348, 0.01, 0.010001, 2.0, 0.01, 1.28552, 0.885225, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.245161, 1.82716, 0.533826, 0.384915, 3.0, 0.242759]
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
single_effs[:, :, 2] .= p[1:12]
single_effs[:, :, 4] .= p[13:24]
single_effs[:, :, [1, 3, 5]] = getODEparams(p[25:end], concs);
single_effs[:, 1, [2, 4]] .= single_effs[:, 1, 1]

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

    x = [1, 2]
    y1 = vcat(mean(efcs[1:2, 1]), mean(efcs[3:6, 1])) # control
    y2 = vcat(mean(efcs[1:2, 2]), mean(efcs[3:6, 2])) # base drug alone
    y3 = vcat(mean(efcs[1:2, 6]), mean(efcs[3:6, 6])) # combination at max of drug B 
    y4 = vcat(mean(single_effs[1:2]), mean(single_effs[3:6])) # single
    scatter(
        x,
        y1,
        color = "blue",
        xlabel = "phase",
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
        alpha=0.5,
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
        y4,
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
    xlims!((0.0, 3.0))
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.03, ymax))
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
    p1 = plot_pCombo(combinEffects[1:6, :, 1], single_effs[1:6, end, 1], "Lap", 1.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "A", 0.15)
    p2 = plot_pCombo(combinEffects[7:12, :, 1], single_effs[7:12, end, 1], "Lap", 0.02, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "B", 0.002)
    p3 = plot_pCombo(combinEffects[1:6, :, 2], single_effs[1:6, end, 3], "Gem", 1.5, "Palbo 50 nM", "Palbo 50 nM + Gem", "progression rates [1/hr]", "C", 0.15)
    p4 = plot_pCombo(combinEffects[7:12, :, 2], single_effs[7:12, end, 3], "Gem", 0.02, "Palbo 50 nM", "Palbo 50 nM + Gem", "death rates [1/hr]", "D", 0.002)
    p5 = plot_pCombo(combinEffects[1:6, :, 3], single_effs[1:6, end, 5], "Palbo", 1.5, "Gem 10 nM", "Gem 10 nM + Palbo", "progression rates [1/hr]", "E", 0.15)
    p6 = plot_pCombo(combinEffects[7:12, :, 3], single_effs[7:12, end, 5], "Palbo", 0.02, "Gem 10 nM", "Gem 10 nM + Palbo", "death rates [1/hr]", "F", 0.002)
    p7 = plot_pCombo(combinEffects[1:6, :, 4], single_effs[1:6, end, 1], "Lap", 1.5, "Gem 10 nM", "Gem 10 nM + Lap", "progression rates [1/hr]", "G", 0.15)
    p8 = plot_pCombo(combinEffects[7:12, :, 4], single_effs[7:12, end, 1], "Lap", 0.02, "Gem 10 nM", "Gem 10 nM + Lap", "death rates [1/hr]", "H", 0.002)
    p9 = plot_pCombo(combinEffects[1:6, :, 5], single_effs[1:6, end, 5], "Palbo", 1.5, "Lpt 100 nM", "Lap 100 nM + Palbo", "progression rates [1/hr]", "I", 0.15)
    p10 = plot_pCombo(combinEffects[7:12, :, 5], single_effs[7:12, end, 5], "Palbo", 0.02, "Lpt 100 nM", "Lap 100 nM + Palbo", "death rates [1/hr]", "J", 0.002)
    p11 = plot_pCombo(combinEffects[1:6, :, 6], single_effs[1:6, end, 3], "Gem", 1.5, "Lpt 100 nM", "Lap 100 nM + Gem", "progression rates [1/hr]", "K", 0.15)
    p12 = plot_pCombo(combinEffects[7:12, :, 6], single_effs[7:12, end, 3], "Gem", 0.02, "Lpt 100 nM", "Lap 100 nM + Gem", "death rates [1/hr]", "L", 0.002)
    p13 = plot_pCombo(combinEffects[1:6, :, 7], single_effs[1:6, end, 3], "Gem", 1.5, "Dox 20 nM", "Dox 20 nM + Gem", "progression rates [1/hr]", "M", 0.15)
    p14 = plot_pCombo(combinEffects[7:12, :, 7], single_effs[7:12, end, 3], "Gem", 0.02, "Dox 20 nM", "Dox 20 nM + Gem", "death rates [1/hr]", "N", 0.002)

    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    figS5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, size = (1800, 900), layout = (2, 7))
    savefig(figS5, "figureS5.svg")
end
