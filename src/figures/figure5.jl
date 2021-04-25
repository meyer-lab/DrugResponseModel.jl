""" This plots time series of combination versus their fit. """

p = [66.8293, 2.4197, 0.00123021, 0.00712543, 6.29059e-5, 0.507736, 0.57301, 1.43432, 0.197454, 0.000746317, 0.0139532, 0.000195288, 0.000198413, 0.0579202, 51.0189, 4.91286, 1.9938, 0.0984497, 0.202416, 0.0898113, 0.462294, 1.96682, 0.701686, 1.72932e-5, 0.000264064, 0.0161778, 1.77773e-5, 0.000346317, 7.81101, 2.99285, 0.176786, 0.517642, 0.171136, 0.630778, 1.88104, 1.55652, 2.50663e-6, 4.3372e-5, 0.00893008, 4.08347e-5, 0.210399, 0.00172866, 36.4208, 1.5254, 0.151112, 0.0177697, 2.56792e-5, 0.000100105, 0.340996, 1.99177, 2.41582e-5, 2.16132e-6, 9.34096e-6, 0.0103405, 0.00016328, 0.204885, 28.1124, 2.64429, 1.15148, 0.282533, 0.46855, 0.260182, 5.08504e-5, 0.945885, 0.000964972, 9.55437e-6, 2.5317e-5, 1.68925e-5, 0.0826429, 5.02612e-5, 83.9744, 0.216313, 0.567271, 0.471734, 0.545467, 0.31523, 0.426661, 1.87149, 1.38939, 0.617567, 0.45418, 1.57687, 0.79395, 0.185619, 0.712989, 1.12255, 0.714262, 0.958377, 0.68417, 0.780506];
concens = zeros(6, 6)
concens[:, 1] = [0, 50, 75, 100, 150, 300]
concens[:, 2] = [0, 50, 55, 60, 67, 80]
concens[:, 3] = [0, 10, 35, 110, 260, 260]
concens[:, 4] = [0, 10, 35, 60, 110, 260]
concens[:, 5] = [0, 20, 25, 30, 37, 50]
concens[:, 6] = [0, 100, 125, 200, 350, 350]

combinEffects = DrugResponseModel.getODEparams(p, concens);

x1 = ["control" "palbo 50 nM" "palbo 50 nM + lap25 nM" "palbo 50 nM + lap 50 nM" "palbo 50 nM + lap 100 nM" "palbo 50 nM + lap 250 nM"]
x2 = ["control" "palbo 50 nM" "palbo 50 nM + gem 5 nM" "palbo 50 nM + gem 10 nM" "palbo 50 nM + gem 17 nM" "palbo 50 nM + gem 30 nM"]
x3 = ["control" "gem 10 nM" "gem 10 nM + palbo 25 nM" "gem 10 nM + palbo 100 nM" "gem 10 nM + palbo 250 nM" "gem 10 nM + palbo 250 nM"]
x4 = ["control" "gem 10 nM" "gem 10 nM + lap 25 nM" "gem 10 nM + lap 50 nM" "gem 10 nM + lap 100 nM" "gem 10 nM + lap 250 nM"]
x5 = ["control" "dox 20 nM" "dox 20 nM + gem 5 nM" "dox 20 nM + gem 10 nM" "dox 20 nM + gem 17 nM" "dox 20 nM + gem 30 nM"]
x6 = ["control" "lap 100 nM" "lap 100 nM + palbo 25 nM" "lap 100 nM + palbo 100 nM" "lap 100 nM + palbo 250 nM" "lap 100 nM + palbo 250 nM"]
GC = JLD.load("GC.jld")["GC"]
GCsim = zeros(2, 193, 6, 6)

t = LinRange(0.0, 96.0, 193)
for k = 1:6 # drug number
    for i = 1:6 # concentration number
        GCsim[1, :, i, k], GCsim[2, :, i, k], _ = predict(combinEffects[:, i, k], combinEffects[:, 1, k], t)
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

function figure5()
    p1 = plot_fig5(GCsim[1, :, :, 1], GC[1, :, :, 1], x1, "Palbo 50 - Lap combo", "G1", "A")
    p2 = plot_fig5(GCsim[2, :, :, 1], GC[2, :, :, 1], x1, "Palbo 50 - Lap combo", "G2", "B")
    p3 = plot_fig5(GCsim[1, :, :, 2], GC[1, :, :, 2], x2, "Palbo 50 - Gem combo", "G1", "C")
    p4 = plot_fig5(GCsim[2, :, :, 2], GC[2, :, :, 2], x2, "Palbo 50 - Gem combo", "G2", "D")
    p5 = plot_fig5(GCsim[1, :, :, 3], GC[1, :, :, 3], x3, "Gem 10 - Palbo combo", "G1", "E")
    p6 = plot_fig5(GCsim[2, :, :, 3], GC[2, :, :, 3], x3, "Gem 10 - Palbo combo", "G2", "F")
    p7 = plot_fig5(GCsim[1, :, :, 4], GC[1, :, :, 4], x4, "Gem 10 - Lap combo", "G1", "G")
    p8 = plot_fig5(GCsim[2, :, :, 4], GC[2, :, :, 4], x4, "Gem 10 - Lap combo", "G2", "H")
    p9 = plot_fig5(GCsim[1, :, :, 5], GC[1, :, :, 5], x5, "Dox 20 - Gem combo", "G1", "I")
    p10 = plot_fig5(GCsim[2, :, :, 5], GC[2, :, :, 5], x5, "Dox 20 - Gem combo", "G2", "J")
    p11 = plot_fig5(GCsim[1, :, :, 6], GC[1, :, :, 6], x6, "Lap 100 - Palbo combo", "G1", "K")
    p12 = plot_fig5(GCsim[2, :, :, 6], GC[2, :, :, 6], x6, "Lap 100 - Palbo combo", "G2", "L")

    fig5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (1400, 1050), layout = (3, 4))
    savefig(fig5, "figure5.svg")
end