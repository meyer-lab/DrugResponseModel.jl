""" plot combined effects. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

p=[22.6116, 0.906908, 0.076665, 0.00321998, 0.140348, 0.623161, 2.47439, 2.46818, 0.00107404, 0.00202817, 6.81464e-5, 0.00949839, 0.00309843, 0.130674, 1.11947, 0.174963, 2.27961, 0.0117939, 0.00275518, 0.162717, 4.50623, 4.9279, 0.00698243, 0.00333901, 2.85654e-5, 0.0241219, 0.00038035, 0.00244354, 20.2521, 0.20926, 0.449486, 0.00905644, 0.376322, 0.540614, 2.42192, 2.4884, 3.00388e-5, 1.50455e-5, 0.000278125, 0.000434414, 0.0119472, 0.184269, 0.293064, 0.784923, 0.248949, 0.378405, 0.582056, 2.03156]
concens = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
combinEffects = DrugResponseModel.my_helper(p)
single_effs = DrugResponseModel.getODEparams(p, concens);

function figureS5()
    p1 = plot_pCombo(combinEffects[1:6, :, 1], single_effs[1:6, end, 1], "Lap", 4.25, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "A", 0.35)
    p2 = plot_pCombo(combinEffects[7:12, :, 1], single_effs[1:6, end, 1], "Lap", 0.3, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "B", 0.06)
    p3 = plot_pCombo(combinEffects[1:6, :, 2], single_effs[1:6, end, 2], "Gem", 4.25, "Palbo 50 nM", "Palbo 50 nM + Gem", "progression rates [1/hr]", "C", 0.35)
    p4 = plot_pCombo(combinEffects[7:12, :, 2], single_effs[1:6, end, 2], "Gem", 0.3, "Palbo 50 nM", "Palbo 50 nM + Gem", "death rates [1/hr]", "D", 0.06)
    p5 = plot_pCombo(combinEffects[1:6, :, 3], single_effs[1:6, end, 3], "Palbo", 4.25, "Gem 10 nM", "Gem 10 nM + Palbo", "progression rates [1/hr]", "E", 0.35)
    p6 = plot_pCombo(combinEffects[7:12, :, 3], single_effs[1:6, end, 3], "Palbo", 0.3, "Gem 10 nM", "Gem 10 nM + Palbo", "death rates [1/hr]", "F", 0.06)
    p7 = plot_pCombo(combinEffects[1:6, :, 4], single_effs[1:6, end, 1], "Lap", 4.25, "Gem 10 nM", "Gem 10 nM + Lap", "progression rates [1/hr]", "G", 0.35)
    p8 = plot_pCombo(combinEffects[7:12, :, 4], single_effs[1:6, end, 1], "Lap", 0.5, "Gem 10 nM", "Gem 10 nM + Lap", "death rates [1/hr]", "H", 0.06)
    p9 = plot_pCombo(combinEffects[1:6, :, 5], single_effs[1:6, end, 1], "Lap", 4.25, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "C", 0.35)
    p10 = plot_pCombo(combinEffects[7:12, :, 5], single_effs[7:12, end, 1], "Lap", 0.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "D", 0.06)

    figS5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, size = (1600, 700), layout = (2, 5))
    savefig(figS5, "figureS5.svg")
end
