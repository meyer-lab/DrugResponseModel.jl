""" plot combined effects. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

p = [0.189371, 1.31225, 0.331636, 3.0, 3.0, 0.274225, 3.32364e-9, 1.48818e-9, 1.5106e-9, 2.04528e-9, 2.3648e-8, 0.0132698, 33.8183, 1.15343, 1.16521e-8, 0.0249639, 0.392687, 0.402105, 2.3868e-8, 0.338073, 0.00245421, 1.36638e-8, 1.62531e-9, 9.78796e-8, 0.0592202, 0.0143732, 0.228744, 2.5, 1.23386, 0.347899, 2.49999, 0.210224]
combinEffects = DrugResponseModel.getODEparams(p, concens);

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
