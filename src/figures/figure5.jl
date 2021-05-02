""" This plots time series of combination versus their fit. """

p =  [22.6116, 0.906908, 0.076665, 0.00321998, 0.140348, 0.623161, 2.47439, 2.46818, 0.00107404, 0.00202817, 6.81464e-5, 0.00949839, 0.00309843, 0.130674, 1.11947, 0.174963, 2.27961, 0.0117939, 0.00275518, 0.162717, 4.50623, 4.9279, 0.00698243, 0.00333901, 2.85654e-5, 0.0241219, 0.00038035, 0.00244354, 20.2521, 0.20926, 0.449486, 0.00905644, 0.376322, 0.540614, 2.42192, 2.4884, 3.00388e-5, 1.50455e-5, 0.000278125, 0.000434414, 0.0119472, 0.184269, 0.293064, 0.784923, 0.248949, 0.378405, 0.582056, 2.03156]

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
function plot_fig5(time, g1, g1data, labels, tite, G, subPlabel)

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
        color = "orange",
        label = "Emax $drugB",
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
    p1 = plot_fig5(t, GCsim[1, :, :, 1], GC[1, :, :, 1], x1, "Palbo 50 + Lap combo", "G1", "A")
    p2 = plot_fig5(t, GCsim[2, :, :, 1], GC[2, :, :, 1], x1, "Palbo 50 + Lap combo", "G2", "B")
    p3 = plot_fig5(t, GCsim[1, :, :, 2], GC[1, :, :, 2], x2, "Palbo 50 + Gem combo", "G1", "C")
    p4 = plot_fig5(t, GCsim[2, :, :, 2], GC[2, :, :, 2], x2, "Palbo 50 + Gem combo", "G2", "D")
    p5 = plot_fig5(t, GCsim[1, :, :, 3], GC[1, :, :, 3], x3, "Gem 10 + Palbo combo", "G1", "E")
    p6 = plot_fig5(t, GCsim[2, :, :, 3], GC[2, :, :, 3], x3, "Gem 10 + Palbo combo", "G2", "F")
    p7 = plot_fig5(t, GCsim[1, :, :, 4], GC[1, :, :, 4], x4, "Gem 10 + Lap combo", "G1", "G")
    p8 = plot_fig5(t, GCsim[2, :, :, 4], GC[2, :, :, 4], x4, "Gem 10 + Lap combo", "G2", "H")
    p9 = plot_fig5(t, GCsim[1, :, :, 5], GC[1, :, :, 6], x6, "Lap 100 + Palbo combo", "G1", "I")
    p10 = plot_fig5(t, GCsim[2, :, :, 5], GC[2, :, :, 6], x6, "Lap 100 + Palbo combo", "G2", "J")

    fig5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, size = (1600, 700), layout = (2, 5))
    savefig(fig5, "figure5.svg")
end
