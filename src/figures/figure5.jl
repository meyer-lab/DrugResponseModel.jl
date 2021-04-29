""" This plots time series of combination versus their fit. """

# p = [0.196747, 0.367277, 2.13154, 2.19888, 0.0774228, 2.49933, 0.000542219, 1.22351e-5, 0.000305728, 0.000354421, 0.000134259, 0.0398953, 36.4052, 1.0259, 0.00799788, 0.000938322, 0.156959, 0.259583, 0.495498, 1.82072, 0.0151445, 0.000100429, 0.000397396, 0.00967118, 0.0183582, 0.00186278, 1.01035, 0.187416, 2.39342, 0.230667, 0.00947831, 0.000180932, 0.490201, 1.80919, 0.0157187, 0.000113442, 1.27606e-5, 0.000116357, 0.0218928, 0.0756593, 19.8996, 0.341808, 2.40293, 0.0992919, 0.252608, 2.4763, 0.246814, 2.47739, 0.0264986, 1.21131e-5, 0.00942984, 0.000439138, 0.00212002, 0.0694528, 0.817749, 0.295525, 0.613096, 0.417977, 0.377789, 0.461066]
p =  [42.287, 1.23679, 0.00092215, 0.000171948, 0.105891, 0.106112, 0.575983, 0.43957, 0.0130425, 1.4023e-5, 0.000612624, 6.37907e-5, 0.00556483, 0.0380514, 1.00236, 1.09617, 0.279774, 0.319337, 2.48502, 0.0866805, 0.139491, 0.524595, 1.7475e-5, 2.01682e-5, 0.0128685, 3.66501e-5, 2.81847e-5, 0.0283298, 6.96782, 0.472736, 0.8431, 0.121317, 0.00609095, 2.49417, 0.24863, 0.507952, 0.00717231, 1.10748e-5, 0.00275721, 0.000171418, 3.52874e-5, 0.0182033, 2.49528, 0.234771, 1.69371, 0.761456, 0.372818, 0.233638]
concs = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
x1 = ["control" "palbo 50 nM" "palbo 50 nM + lap25 nM" "palbo 50 nM + lap 50 nM" "palbo 50 nM + lap 100 nM" "palbo 50 nM + lap 250 nM"]
x2 = ["control" "palbo 50 nM" "palbo 50 nM + gem 5 nM" "palbo 50 nM + gem 10 nM" "palbo 50 nM + gem 17 nM" "palbo 50 nM + gem 30 nM"]
x3 = ["control" "gem 10 nM" "gem 10 nM + palbo 25 nM" "gem 10 nM + palbo 50 nM" "gem 10 nM + palbo 100 nM" "gem 10 nM + palbo 250 nM"]
x4 = ["control" "gem 10 nM" "gem 10 nM + lap 25 nM" "gem 10 nM + lap 50 nM" "gem 10 nM + lap 100 nM" "gem 10 nM + lap 250 nM"]
x5 = ["control" "dox 20 nM" "dox 20 nM + gem 5 nM" "dox 20 nM + gem 10 nM" "dox 20 nM + gem 17 nM" "dox 20 nM + gem 30 nM"]
x6 = ["control" "lap 100 nM" "lap 100 nM + palbo 25 nM" "lap 100 nM + palbo 50 nM" "lap 100 nM + palbo 100 nM" "lap 100 nM + palbo 250 nM"]
GC = JLD.load("GC.jld")["GC"]
GCsim = zeros(2, 193, 6, 5)

t = LinRange(0.0, 96.0, 193)

combinEffects = DrugResponseModel.my_helper(p)
for j=1:5
    for i = 1:6 # concentration number
        GCsim[1, :, i, j], GCsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
    end
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

function plot_pCombo(efcs, single_effs, ymax, label2, Phasename, ylabel, subPlabel, plus)

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
    scatter!(
        x,
        single_effs,
        color = "orange",
        label = "Emax Lap",
        alpha=0.6,
        markerstrokewidth = 0,
        markershape=:xcross,
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
    p1 = plot_fig5(GCsim[1, :, :, 1], GC[1, :, :, 1], x1, "Palbo 50 + Lap combo", "G1", "A")
    p2 = plot_fig5(GCsim[2, :, :, 1], GC[2, :, :, 1], x1, "Palbo 50 + Lap combo", "G2", "B")
    p3 = plot_fig5(GCsim[1, :, :, 2], GC[1, :, :, 2], x2, "Palbo 50 + Gem combo", "G1", "C")
    p4 = plot_fig5(GCsim[2, :, :, 2], GC[2, :, :, 2], x2, "Palbo 50 + Gem combo", "G2", "D")
    p5 = plot_fig5(GCsim[1, :, :, 3], GC[1, :, :, 3], x3, "Gem 10 + Palbo combo", "G1", "E")
    p6 = plot_fig5(GCsim[2, :, :, 3], GC[2, :, :, 3], x3, "Gem 10 + Palbo combo", "G2", "F")
    p7 = plot_fig5(GCsim[1, :, :, 4], GC[1, :, :, 4], x4, "Gem 10 + Lap combo", "G1", "G")
    p8 = plot_fig5(GCsim[2, :, :, 4], GC[2, :, :, 4], x4, "Gem 10 + Lap combo", "G2", "H")
    # p9 = plot_fig5(GCsim[1, :, :, 5], GC[1, :, :, 5], x5, "Dox 20 + Gem combo", "G1", "I")
    # p10 = plot_fig5(GCsim[2, :, :, 5], GC[2, :, :, 5], x5, "Dox 20 + Gem combo", "G2", "J")
    p9 = plot_fig5(GCsim[1, :, :, 5], GC[1, :, :, 6], x6, "Lap 100 + Palbo combo", "G1", "I")
    p10 = plot_fig5(GCsim[2, :, :, 5], GC[2, :, :, 6], x6, "Lap 100 + Palbo combo", "G2", "J")
    # p3 = plot_pCombo(combinEffects[1:6, :], single_effs[1:6, end, 1], 3.75, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "C", 0.35)
    # p4 = plot_pCombo(combinEffects[7:12, :], single_effs[7:12, end, 1], 0.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "D", 0.06)

    fig5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, size = (1650, 850), layout = (2, 5))
    savefig(fig5, "figure51.svg")
end
