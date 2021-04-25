""" plot combined effects. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

concens = zeros(6, 6)
concens[:, 1] = [0, 50, 75, 100, 150, 300]
concens[:, 2] = [0, 50, 55, 60, 67, 80]
concens[:, 3] = [0, 10, 35, 110, 260, 260]
concens[:, 4] = [0, 10, 35, 60, 110, 260]
concens[:, 5] = [0, 20, 25, 30, 37, 50]
concens[:, 6] = [0, 100, 125, 200, 350, 350]

p = [66.8293, 2.4197, 0.00123021, 0.00712543, 6.29059e-5, 0.507736, 0.57301, 1.43432, 0.197454, 0.000746317, 0.0139532, 0.000195288, 0.000198413, 0.0579202, 51.0189, 4.91286, 1.9938, 0.0984497, 0.202416, 0.0898113, 0.462294, 1.96682, 0.701686, 1.72932e-5, 0.000264064, 0.0161778, 1.77773e-5, 0.000346317, 7.81101, 2.99285, 0.176786, 0.517642, 0.171136, 0.630778, 1.88104, 1.55652, 2.50663e-6, 4.3372e-5, 0.00893008, 4.08347e-5, 0.210399, 0.00172866, 36.4208, 1.5254, 0.151112, 0.0177697, 2.56792e-5, 0.000100105, 0.340996, 1.99177, 2.41582e-5, 2.16132e-6, 9.34096e-6, 0.0103405, 0.00016328, 0.204885, 28.1124, 2.64429, 1.15148, 0.282533, 0.46855, 0.260182, 5.08504e-5, 0.945885, 0.000964972, 9.55437e-6, 2.5317e-5, 1.68925e-5, 0.0826429, 5.02612e-5, 83.9744, 0.216313, 0.567271, 0.471734, 0.545467, 0.31523, 0.426661, 1.87149, 1.38939, 0.617567, 0.45418, 1.57687, 0.79395, 0.185619, 0.712989, 1.12255, 0.714262, 0.958377, 0.68417, 0.780506];

combinEffects = DrugResponseModel.getODEparams(p, concens);

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
        y2,
        color = "red",
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

function figureS5()
    p1 = plot_pCombo(combinEffects[1:6, :, 1], 2.25, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "A", 0.35)
    p2 = plot_pCombo(combinEffects[7:12, :, 1], 0.2, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "B", 0.06)
    p3 = plot_pCombo(combinEffects[1:6, :, 2], 2.25, "Palbo 50 nM", "Palbo 50 nM + Gem", "progression rates [1/hr]", "C", 0.35)
    p4 = plot_pCombo(combinEffects[7:12, :, 2], 0.2, "Palbo 50 nM", "Palbo 50 nM + Gem", "death rates [1/hr]", "D", 0.06)
    p5 = plot_pCombo(combinEffects[1:6, :, 3], 2.25, "Gem 10 nM", "Gem 10 nM + Palbo", "progression rates [1/hr]", "E", 0.35)
    p6 = plot_pCombo(combinEffects[7:12, :, 3], 0.2, "Gem 10 nM", "Gem 10 nM + Palbo", "death rates [1/hr]", "F", 0.06)
    p7 = plot_pCombo(combinEffects[1:6, :, 4], 2.25, "Gem 10 nM", "Gem 10 nM + Lap", "progression rates [1/hr]", "G", 0.35)
    p8 = plot_pCombo(combinEffects[7:12, :, 4], 0.2, "Gem 10 nM", "Gem 10 nM + Lap", "death rates [1/hr]", "H", 0.06)
    p9 = plot_pCombo(combinEffects[1:6, :, 5], 2.25, "Dox 20 nM", "Dox 20 nM + Gem", "progression rates [1/hr]", "I", 0.35)
    p10 = plot_pCombo(combinEffects[7:12, :, 5], 0.2, "Dox 20 nM", "Dox 20 nM + Gem", "death rates [1/hr]", "J", 0.06)
    p11 = plot_pCombo(combinEffects[1:6, :, 6], 2.25, "Lap 100 nM", "Lap 100 nM + Palbo", "progression rates [1/hr]", "K", 0.35)
    p12 = plot_pCombo(combinEffects[7:12, :, 6], 0.2, "Lap 100 nM", "Lap 100 nM + Palbo", "death rates [1/hr]", "L", 0.06)
    figS5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (1400, 1050), layout = (3, 4))
    savefig(fig5, "figureS5.svg")
end
