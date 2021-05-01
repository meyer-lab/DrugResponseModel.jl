""" plot combined effects. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

p =  [42.287, 1.23679, 0.00092215, 0.000171948, 0.105891, 0.106112, 0.575983, 0.43957, 0.0130425, 1.4023e-5, 0.000612624, 6.37907e-5, 0.00556483, 0.0380514, 1.00236, 1.09617, 0.279774, 0.319337, 2.48502, 0.0866805, 0.139491, 0.524595, 1.7475e-5, 2.01682e-5, 0.0128685, 3.66501e-5, 2.81847e-5, 0.0283298, 6.96782, 0.472736, 0.8431, 0.121317, 0.00609095, 2.49417, 0.24863, 0.507952, 0.00717231, 1.10748e-5, 0.00275721, 0.000171418, 3.52874e-5, 0.0182033, 2.49528, 0.234771, 1.69371, 0.761456, 0.372818, 0.233638]
concens = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
combinEffects = DrugResponseModel.my_helper(p)
single_effs = DrugResponseModel.getODEparams(p, concens);

function figureS5()
    p1 = plot_pCombo(combinEffects[1:6, :, 1], single_effs[1:6, end, 1], "Lap", 3.75, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "A", 0.35)
    p2 = plot_pCombo(combinEffects[7:12, :, 1], single_effs[1:6, end, 1], "Lap", 0.2, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "B", 0.06)
    p3 = plot_pCombo(combinEffects[1:6, :, 2], single_effs[1:6, end, 2], "Gem", 3.75, "Palbo 50 nM", "Palbo 50 nM + Gem", "progression rates [1/hr]", "C", 0.35)
    p4 = plot_pCombo(combinEffects[7:12, :, 2], single_effs[1:6, end, 2], "Gem", 0.2, "Palbo 50 nM", "Palbo 50 nM + Gem", "death rates [1/hr]", "D", 0.06)
    p5 = plot_pCombo(combinEffects[1:6, :, 3], single_effs[1:6, end, 3], "Palbo", 3.75, "Gem 10 nM", "Gem 10 nM + Palbo", "progression rates [1/hr]", "E", 0.35)
    p6 = plot_pCombo(combinEffects[7:12, :, 3], single_effs[1:6, end, 3], "Palbo", 0.2, "Gem 10 nM", "Gem 10 nM + Palbo", "death rates [1/hr]", "F", 0.06)
    p7 = plot_pCombo(combinEffects[1:6, :, 4], single_effs[1:6, end, 1], "Lap", 3.75, "Gem 10 nM", "Gem 10 nM + Lap", "progression rates [1/hr]", "G", 0.35)
    p8 = plot_pCombo(combinEffects[7:12, :, 4], single_effs[1:6, end, 1], "Lap", 0.5, "Gem 10 nM", "Gem 10 nM + Lap", "death rates [1/hr]", "H", 0.06)
    p9 = plot_pCombo(combinEffects[1:6, :, 5], single_effs[1:6, end, 1], "Lap", 3.75, "Palbo 50 nM", "Palbo 50 nM + Lap", "progression rates [1/hr]", "C", 0.35)
    p10 = plot_pCombo(combinEffects[7:12, :, 5], single_effs[7:12, end, 1], "Lap", 0.5, "Palbo 50 nM", "Palbo 50 nM + Lap", "death rates [1/hr]", "D", 0.06)

    figS5 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, size = (1600, 700), layout = (2, 5))
    savefig(figS5, "figureS5.svg")
end
