""" Drug combination for new data. """

function combination_newdrugs()

    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    # bliss on end-point cell numbers
    total = tensor[1, :, :, :] .+ tensor[2, :, :, :]
    # find combination of total cell numbers at the end time point
    combin = DrugResponseModel.blissCellNum(total[182, :, :]) # 8 x 8 x 15

    # bliss from model predictions
    params = [9.20945, 2.12023, 2.98291, 0.430111, 0.09537, 0.0956342, 0.242984, 0.0828915, 1.71385, 0.299598, 0.000505937, 0.0647114, 5.01559e-6, 1.49885e-5, 0.00544226, 4.35115e-6, 3.03577e-5, 2.15171e-5, 13.9863, 2.77415, 3.06754, 3.09629, 0.129791, 0.002272, 2.12184, 0.461891, 0.362844, 0.0293403, 0.000990191, 0.00063532, 2.98906e-5, 3.40495e-6, 0.000103298, 5.23445e-7, 7.58079e-6, 0.00361932, 1.11528, 1.93197, 3.82999, 3.8369, 7.05901e-6, 3.52973e-5, 0.714193, 0.0650093, 0.521932, 0.0465728, 0.00111718, 6.65118e-6, 0.00254141, 0.0378996, 9.51571e-5, 7.83845e-7, 0.000242097, 0.00102866, 1.5658, 3.63417, 3.88274, 3.97669, 0.0956175, 0.64807, 2.13597, 0.396049, 1.0315, 0.395203, 0.00162117, 0.00053917, 0.000125216, 0.0929912, 0.00102563, 0.000176756, 0.00184695, 0.0261891, 13.4388, 0.395591, 0.0829676, 0.000898367, 0.0907073, 0.000483117, 0.0306672, 0.166496, 0.000362649, 0.399909, 0.00273783, 0.0227258, 3.30817e-6, 8.66783e-6, 9.05333e-6, 2.4177e-5, 0.000184217, 2.80325e-6, 3.56209, 1.07825, 3.95194, 0.447834, 0.0389113, 0.110297, 2.13567, 0.27307, 3.36894, 0.368404, 0.000296116, 0.0415469, 4.95713e-6, 2.47911e-5, 2.06936e-5, 3.39401e-5, 0.0908003, 3.87974e-6, 9.97915, 9.99163, 0.288515, 0.632015, 2.12838, 0.320644, 6.53941, 0.325873]
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
                    g1, g2, _ = DrugResponseModel.predict(param_bliss[:, h, u], param_bliss[:, 1, 1], 92.0)
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