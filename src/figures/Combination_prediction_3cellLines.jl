""" Drug combination for new data. """

function combination_newdrugs()

    tensor, names, concs, conds = DrugResponseModel.hcc_all()
    # bliss on end-point cell numbers
    total = tensor[1, :, :, :] .+ tensor[2, :, :, :]
    # find combination of total cell numbers at the end time point
    combin = DrugResponseModel.blissCellNum(total[182, :, :]) # 8 x 8 x 15

    # bliss from model predictions
    # this is from HCC1143. Please find parameters from other two cell lines from src/figures/figureS5.jl
    params = [34.9432, 27.829, 0.291144, 0.0244977, 2.64026, 0.115686, 0.442475, 3.61059e-5, 0.422046, 0.330667, 9.47844e-6, 3.39363e-6, 0.000705232, 2.34289e-5, 1.23165e-5, 2.20203e-6, 2.92093, 8.09313e-6, 159.166, 13.0877, 1.5532, 0.00387504, 3.63611, 0.260283, 0.450801, 0.0250362, 0.367939, 0.252831, 0.00107609, 0.000115707, 0.00260766, 6.82032e-5, 4.10667e-5, 9.84008e-6, 1.46092, 0.000134889, 0.741738, 2.4431, 3.9629, 0.000761489, 3.94527, 0.348368, 0.331174, 2.27071, 0.0737219, 0.591972, 0.00185578, 0.0162289, 0.00154473, 2.32892e-5, 1.1137e-7, 5.03267e-5, 7.37255e-7, 0.00531612, 1.46012, 4.8298, 0.517965, 0.08571, 0.736769, 0.0322651, 0.290752, 3.31035, 0.0641698, 0.292762, 9.08443e-5, 5.09103e-5, 0.000280376, 0.00737825, 0.000304891, 0.00129649, 3.0795e-5, 4.0636e-5, 86.4606, 0.338272, 3.88389, 0.00102644, 3.58318, 0.000195763, 0.104922, 3.96811, 0.24616, 0.657727, 0.000600609, 0.0132196, 0.000113768, 3.56294e-6, 4.42152e-6, 4.73675e-5, 8.75234e-6, 2.35412e-5, 3.716, 1.19482, 0.564743, 0.0361228, 3.98035, 0.110105, 0.403263, 3.96816, 0.421521, 0.640074, 0.0878646, 0.000286413, 0.000384056, 1.16993e-5, 1.07933e-5, 0.00215579, 8.21145e-5, 0.000489024, 3.98591, 0.588244, 3.98874, 0.268013, 0.441755, 3.99003, 0.365446, 0.505571]
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end
    pode = getODEparams(params, cs)
    model_combin = zeros(8, 8, 15)
    x = 1
    num = 1
    for i = 1:5
        for k = (i + 1):6
            param_bliss = DrugResponseModel.AllBliss_params(pode[:, :, i], pode[:, :, k])
            for h = 1:8
                for u = 1:8
                    g1, g2, _ = DrugResponseModel.predict(param_bliss[:, h, u], param_bliss[:, 1, 1], 190.0)
                    model_combin[h, u, x] = g1 .+ g2
                end
            end
            x += 1
        end
        num += 1
    end

    return combin, model_combin, cs
end

function for_each_pair(conc_1, conc_2, combin, model_combin, dr1, dr2)
    plt = []
    for i=1:8
        df = [DataFrame(x=log10.(conc_1), y=combin[:, i], label="on cell num"), DataFrame(x=log10.(conc_1), y=model_combin[:, i], label="on rates")]
        DF = vcat(df...)

        c2 = conc_2[i]
        push!(plt, Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.title(string(dr1, " + $c2 nM ", dr2)), Guide.xlabel(string("concentration [nM] ", dr1)), style(line_width=1mm), Coord.Cartesian(ymin=0.0,ymax=4.0)))
    end
    return plt
end

function figure10()
    setGadflyTheme()
    combin, model_combin, cs = combination_newdrugs()
    p1 = for_each_pair(cs[:, 1], cs[:, 2], combin[:, :, 1], model_combin[:, :, 1], "BEZ", "DOX")
    p2 = for_each_pair(cs[:, 1], cs[:, 3], combin[:, :, 2], model_combin[:, :, 2], "BEZ", "GEM")
    p3 = for_each_pair(cs[:, 1], cs[:, 4], combin[:, :, 3], model_combin[:, :, 3], "BEZ", "Taxol")
    p4 = for_each_pair(cs[:, 1], cs[:, 5], combin[:, :, 4], model_combin[:, :, 4], "BEZ", "Palbo")
    p5 = for_each_pair(cs[:, 1], cs[:, 6], combin[:, :, 5], model_combin[:, :, 5], "BEZ", "Trametinib")

    pl1 = plotGrid((5, 8), [p1..., p2..., p3..., p4..., p5...];)
    return draw(SVG("figure10_bez.svg", 25inch, 18inch), pl1)
end

function figure101()
    setGadflyTheme()
    combin, model_combin, cs = combination_newdrugs()
    p1 = for_each_pair(cs[:, 2], cs[:, 3], combin[:, :, 6], model_combin[:, :, 6], "DOX", "GEM")
    p2 = for_each_pair(cs[:, 2], cs[:, 4], combin[:, :, 7], model_combin[:, :, 7], "DOX", "Taxol")
    p3 = for_each_pair(cs[:, 2], cs[:, 5], combin[:, :, 8], model_combin[:, :, 8], "DOX", "Palbo")
    p4 = for_each_pair(cs[:, 2], cs[:, 6], combin[:, :, 9], model_combin[:, :, 9], "DOX", "Trametinib")

    pl1 = plotGrid((4, 8), [p1..., p2..., p3..., p4...];)
    return draw(SVG("figure10_dox.svg", 25inch, 12inch), pl1)
end

function figure102()
    setGadflyTheme()
    combin, model_combin, cs = combination_newdrugs()
    p1 = for_each_pair(cs[:, 3], cs[:, 4], combin[:, :, 10], model_combin[:, :, 10], "GEM", "Taxol")
    p2 = for_each_pair(cs[:, 3], cs[:, 5], combin[:, :, 11], model_combin[:, :, 11], "GEM", "Palbo")
    p3 = for_each_pair(cs[:, 3], cs[:, 6], combin[:, :, 12], model_combin[:, :, 12], "GEM", "Trametinib")

    pl1 = plotGrid((3, 8), [p1..., p2..., p3...];)
    return draw(SVG("figure10_gem.svg", 25inch, 8inch), pl1)
end

function figure103()
    setGadflyTheme()
    combin, model_combin, cs = combination_newdrugs()
    p1 = for_each_pair(cs[:, 4], cs[:, 5], combin[:, :, 13], model_combin[:, :, 13], "Taxol", "Palbo")
    p2 = for_each_pair(cs[:, 4], cs[:, 6], combin[:, :, 14], model_combin[:, :, 14], "Taxol", "Trametinib")
    p3 = for_each_pair(cs[:, 5], cs[:, 6], combin[:, :, 15], model_combin[:, :, 15], "Palbo", "Trametinib")

    pl1 = plotGrid((3, 8), [p1..., p2..., p3...];)
    return draw(SVG("figure10_taxol_palbo.svg", 25inch, 8inch), pl1)
end