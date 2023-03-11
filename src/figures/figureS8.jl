""" Supplementary Figure 8: the quantified effects of drugs for the three cell lines (HCC1143, MDAMB157, 21MT1) when fitting the drugs at once. """

""" Helper function for plotting the effects. """
function effectsplot_helper(params, cs)
    efcs = getODEparams(params, cs)
    # params at EC50
    ec50 = zeros(16, 6)
    conc_ec50 = zeros((1, 6))
    conc_ec50[1, 1:6] = [params[18*k+1] for k = 0:5]
    ec50 = getODEparams(params, conc_ec50)[:, 1, :]
    print(conc_ec50)

    # phase durations
    # @ control
    gi = zeros(2, 2, 6)
    gi[1, 1, :] .= (2 ./ efcs[1, 1, :] .+ 2 ./ efcs[2, 1, :] .+ 2 ./ efcs[3, 1, :] .+ 2 ./ efcs[4, 1, :])
    gi[2, 1, :] .= (5 ./ efcs[5, 1, :] .+ 5 ./ efcs[6, 1, :] .+ 5 ./ efcs[7, 1, :] .+ 5 ./ efcs[8, 1, :])

    # @ ec50
    gi[1, 2, :] .= (2 ./ ec50[1, :] .+ 2 ./ ec50[2, :] .+ 2 ./ ec50[3, :] .+ 2 ./ ec50[4, :])
    gi[2, 2, :] .= (5 ./ ec50[5, :] .+ 5 ./ ec50[6, :] .+ 5 ./ ec50[7, :] .+ 5 ./ ec50[8, :])

    # cell deaths
    deathContG1 = zeros(4, 6)
    deathEC50G1 = zeros(4, 6)
    deathContG1[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[1, 1, :])
    deathEC50G1[1, :] .= (ec50[9, :]) ./ (ec50[9, :] .+ ec50[1, :])
    for i = 2:4
        deathContG1[i, :] .= (1 .- deathContG1[i - 1, :]) .* (efcs[i + 8, 1, :]) ./ (efcs[i + 8, 1, :] .+ efcs[i, 1, :])
        deathEC50G1[i, :] .= (1 .- deathEC50G1[i - 1, :]) .* (ec50[i + 8, :]) ./ (ec50[i + 8, :] .+ ec50[i, :])
    end
    deathContG2 = zeros(4, 6)
    deathEC50G2 = zeros(4, 6)
    deathContG2[1, :] .= (efcs[13, 1, :]) ./ (efcs[13, 1, :] .+ efcs[5, 1, :])
    deathEC50G2[1, :] .= (ec50[13, :]) ./ (ec50[13, :] .+ ec50[5, :])
    for i = 14:16
        deathContG2[i - 12, :] = (1 .- deathContG2[i - 13, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 8, 1, :])
        deathEC50G2[i - 12, :] = (1 .- deathEC50G2[i - 13, :]) .* (ec50[i, :]) ./ (ec50[i, :] .+ ec50[i - 8, :])
    end

    return gi, deathContG1, deathEC50G1, deathContG2, deathEC50G2

end


""" Plotting the effects of all three cell lines on top of each other on a plot. """
function figure800()

    ENV["GKSwstype"]="nul"
    _, _, hcc_concs, _ = DrugResponseModel.hcc_all()
    _, _, mt1_concs, _ = DrugResponseModel.mt1_all()
    _, _, mdamb_concs, _ = DrugResponseModel.mda_all()

    num_drugs = length(hcc_concs)
    cs = zeros(8, num_drugs, 3)
    for i=1:num_drugs
        cs[:, i, 1] .= hcc_concs[i]
        cs[:, i, 2] .= mt1_concs[i]
        cs[:, i, 3] .= mdamb_concs[i]
    end

    # 6-drug fits:
    params_hcc = [34.9432, 27.829, 0.291144, 0.0244977, 2.64026, 0.115686, 0.442475, 3.61059e-5, 0.422046, 0.330667, 9.47844e-6, 3.39363e-6, 0.000705232, 2.34289e-5, 1.23165e-5, 2.20203e-6, 2.92093, 8.09313e-6, 159.166, 13.0877, 1.5532, 0.00387504, 3.63611, 0.260283, 0.450801, 0.0250362, 0.367939, 0.252831, 0.00107609, 0.000115707, 0.00260766, 6.82032e-5, 4.10667e-5, 9.84008e-6, 1.46092, 0.000134889, 0.741738, 2.4431, 3.9629, 0.000761489, 3.94527, 0.348368, 0.331174, 2.27071, 0.0737219, 0.591972, 0.00185578, 0.0162289, 0.00154473, 2.32892e-5, 1.1137e-7, 5.03267e-5, 7.37255e-7, 0.00531612, 1.46012, 4.8298, 0.517965, 0.08571, 0.736769, 0.0322651, 0.290752, 3.31035, 0.0641698, 0.292762, 9.08443e-5, 5.09103e-5, 0.000280376, 0.00737825, 0.000304891, 0.00129649, 3.0795e-5, 4.0636e-5, 86.4606, 0.338272, 3.88389, 0.00102644, 3.58318, 0.000195763, 0.104922, 3.96811, 0.24616, 0.657727, 0.000600609, 0.0132196, 0.000113768, 3.56294e-6, 4.42152e-6, 4.73675e-5, 8.75234e-6, 2.35412e-5, 3.716, 1.19482, 0.564743, 0.0361228, 3.98035, 0.110105, 0.403263, 3.96816, 0.421521, 0.640074, 0.0878646, 0.000286413, 0.000384056, 1.16993e-5, 1.07933e-5, 0.00215579, 8.21145e-5, 0.000489024, 3.98591, 0.588244, 3.98874, 0.268013, 0.441755, 3.99003, 0.365446, 0.505571]
    params_21mt1 = [1.55185, 1.45433, 3.99842, 0.00724615, 2.46356, 0.102305, 2.60303, 3.3045, 1.04112, 2.71372, 1.82752e-5, 0.0057077, 1.04167e-5, 0.445341, 5.64574e-5, 3.24302e-5, 8.55308e-7, 3.00144e-5, 19.1045, 0.698229, 1.04623, 0.086005, 1.00304, 1.84451e-6, 2.55208, 3.29383, 0.606597, 2.67352, 3.3935e-6, 5.37281e-7, 0.000554417, 0.00117723, 1.21694e-5, 7.34184e-5, 1.15653e-5, 4.54135e-6, 220.841, 14.6972, 3.83869, 0.0975135, 2.15286, 0.302349, 2.5924, 3.29021, 0.438998, 2.65367, 0.000445606, 7.23863e-6, 0.00012159, 1.74178e-5, 1.34838e-5, 5.1107e-5, 3.84569e-6, 3.03796e-5, 12.5217, 2.23514, 0.161317, 0.249599, 2.47102, 0.0467338, 0.0740153, 1.24226, 0.495188, 0.67336, 3.46356e-7, 3.85531e-6, 4.22593e-5, 9.41142e-8, 6.0454e-7, 2.56087e-7, 4.66992e-7, 2.07443e-7, 51.6852, 2.10108, 0.0332774, 0.00824266, 0.00151619, 0.000457234, 2.60334, 3.28454, 0.0951307, 2.10065e-5, 0.0132314, 0.000176011, 0.646514, 3.2941, 9.48263e-5, 1.80983e-6, 1.39552e-6, 2.5178e-7, 1.51555, 2.32919, 1.65229, 0.0336268, 0.329872, 0.000369634, 2.53846, 4.65599e-5, 0.420307, 4.96256e-5, 4.27684e-6, 2.02799e-6, 1.51012e-7, 5.99119e-7, 2.30595e-5, 3.27434e-7, 3.62375e-6, 9.15975e-8, 3.99449, 0.603004, 2.4585, 3.99946, 2.5756, 3.28379, 0.512648, 2.64776]
    params_mdamb = [2.01648, 2.7393, 0.381298, 0.282195, 0.0426087, 0.00438904, 3.99969, 0.38872, 3.99991, 2.56698e-7, 2.74086e-7, 3.53968e-7, 4.66013e-7, 0.31455, 8.04555e-7, 1.79967e-7, 1.42038e-8, 0.0394145, 3594.01, 38.8183, 3.97541, 3.98819, 0.0427238, 1.01873, 0.0825609, 0.365346, 3.99997, 0.237293, 9.11982e-8, 2.1046e-7, 1.09047e-8, 1.41151e-8, 6.54729e-9, 2.89062e-8, 1.86293e-6, 1.22863e-8, 658.205, 2.11387, 3.99996, 0.150384, 0.0246126, 4.0, 3.7199, 0.296986, 3.99558, 0.619582, 9.13771e-8, 4.40758e-8, 2.34502e-8, 0.804942, 3.73842e-7, 1.71125e-8, 8.63841e-7, 3.34926e-8, 10.9042, 4.84963, 3.99996, 0.0861652, 0.0428161, 0.24123, 0.0676982, 0.127572, 3.1579, 0.367814, 4.89605e-7, 0.000408088, 1.35047e-8, 8.0359e-7, 1.64854e-7, 7.10794e-9, 1.41167e-7, 2.52954e-8, 247.586, 42.2574, 0.699003, 0.0740976, 0.00144902, 3.9916, 0.0833487, 0.338742, 3.97948, 0.175912, 7.82846e-7, 7.6491e-8, 0.0300662, 9.54405e-8, 2.03674e-6, 9.53742e-8, 1.03249e-6, 3.34538e-8, 0.642888, 7.62552, 3.99996, 0.117966, 0.344893, 0.236656, 0.080841, 0.306802, 3.9999, 0.807695, 3.52664e-7, 1.25166e-8, 7.37586e-8, 4.40688e-9, 1.28515e-8, 0.0340858, 1.90609e-7, 8.55822e-9, 3.99998, 3.99998, 0.123506, 3.99037, 3.99998, 0.270454, 4.0, 0.780807]
    gi_h, deathContG1_h, deathEC50G1_h, deathContG2_h, deathEC50G2_h = effectsplot_helper(params_hcc, cs[:, :, 1])
    gi_t, deathContG1_t, deathEC50G1_t, deathContG2_t, deathEC50G2_t = effectsplot_helper(params_21mt1, cs[:, :, 2])
    gi_m, deathContG1_m, deathEC50G1_m, deathContG2_m, deathEC50G2_m = effectsplot_helper(params_mdamb, cs[:, :, 3])

    # x = ["Paclitaxel", "Palbociclib", "Doxorubicin", "Gemcitabine"]
    x = ["Paclitaxel", "Palbociclib", "Trametinib", "BEZ235", "Doxorubicin", "Gemcitabine"]
    dfs1 = DataFrame(x=repeat(x, 3), y=vcat(gi_h[1, 1, :], gi_t[1, 1, :], gi_m[1, 1, :]), y2=vcat(gi_h[1, 2, :], gi_t[1, 2, :], gi_m[1, 2, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=6), label="G1")
    dfs2 = DataFrame(x=repeat(x, 3), y=vcat(gi_h[2, 1, :], gi_t[2, 1, :], gi_m[2, 1, :]), y2=vcat(gi_h[2, 2, :], gi_t[2, 2, :], gi_m[2, 2, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=6), label="SG2")

    dfs3 = DataFrame(x=repeat(x, 3), y=vcat(sum(deathContG1_h, dims = 1)[1, :], sum(deathContG1_t, dims = 1)[1, :], sum(deathContG1_m, dims = 1)[1, :]), y2=vcat(sum(deathEC50G1_h, dims = 1)[1, :], sum(deathEC50G1_t, dims = 1)[1, :], sum(deathEC50G1_m, dims = 1)[1, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=6), label="G1")
    dfs4 = DataFrame(x=repeat(x, 3), y=vcat(sum(deathContG2_h, dims = 1)[1, :], sum(deathContG2_t, dims = 1)[1, :], sum(deathContG2_m, dims = 1)[1, :]), y2=vcat(sum(deathEC50G2_h, dims = 1)[1, :], sum(deathEC50G2_t, dims = 1)[1, :], sum(deathEC50G2_m, dims = 1)[1, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=6), label="SG2")

    p = [Gadfly.plot(dfs1, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),

    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]")),

    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),

    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"))]

    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("SupplementaryFigure8.svg", 10inch, 6inch), pl1)
end
