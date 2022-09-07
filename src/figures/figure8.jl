""" Figure 8: the quantified effects of drugs for HCC cell line when fitting at once. """
function plt_2(i, j, eff_name, ymax)
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end

    params = [9.20945, 2.12023, 2.98291, 0.430111, 0.09537, 0.0956342, 0.242984, 0.0828915, 1.71385, 0.299598, 0.000505937, 0.0647114, 5.01559e-6, 1.49885e-5, 0.00544226, 4.35115e-6, 3.03577e-5, 2.15171e-5, 13.9863, 2.77415, 3.06754, 3.09629, 0.129791, 0.002272, 2.12184, 0.461891, 0.362844, 0.0293403, 0.000990191, 0.00063532, 2.98906e-5, 3.40495e-6, 0.000103298, 5.23445e-7, 7.58079e-6, 0.00361932, 1.11528, 1.93197, 3.82999, 3.8369, 7.05901e-6, 3.52973e-5, 0.714193, 0.0650093, 0.521932, 0.0465728, 0.00111718, 6.65118e-6, 0.00254141, 0.0378996, 9.51571e-5, 7.83845e-7, 0.000242097, 0.00102866, 1.5658, 3.63417, 3.88274, 3.97669, 0.0956175, 0.64807, 2.13597, 0.396049, 1.0315, 0.395203, 0.00162117, 0.00053917, 0.000125216, 0.0929912, 0.00102563, 0.000176756, 0.00184695, 0.0261891, 13.4388, 0.395591, 0.0829676, 0.000898367, 0.0907073, 0.000483117, 0.0306672, 0.166496, 0.000362649, 0.399909, 0.00273783, 0.0227258, 3.30817e-6, 8.66783e-6, 9.05333e-6, 2.4177e-5, 0.000184217, 2.80325e-6, 3.56209, 1.07825, 3.95194, 0.447834, 0.0389113, 0.110297, 2.13567, 0.27307, 3.36894, 0.368404, 0.000296116, 0.0415469, 4.95713e-6, 2.47911e-5, 2.06936e-5, 3.39401e-5, 0.0908003, 3.87974e-6, 9.97915, 9.99163, 0.288515, 0.632015, 2.12838, 0.320644, 6.53941, 0.325873]
    pp = getODEparams(params, cs)

    dfs = [DataFrame(x=concs[i], y=pp[Int(j+k-1), :, i], label=string(eff_name*"$k")) for k = 1:4]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.title(names[i]), Guide.xlabel("concentrations [nM]"), Guide.ylabel(eff_name), style(line_width=1mm), Coord.Cartesian(ymin=0.0,ymax=ymax))
end

function figure8()
    setGadflyTheme()
    plts = []
    for i=1:6
        push!(plts, plt_2(i, 1, "G1 progression", 4.0))
        push!(plts, plt_2(i, 5, "S-G2 progression", 4.0))
        push!(plts, plt_2(i, 9, "G1 cell death", 0.5))
        push!(plts, plt_2(i, 12, "S-G2 cell death", 0.5))
    end

    pl1 = plotGrid((6, 4), [plts...];)
    return draw(SVG("figure8.svg", 16inch, 20inch), pl1)
end

function figure81()
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end

    params = [9.20945, 2.12023, 2.98291, 0.430111, 0.09537, 0.0956342, 0.242984, 0.0828915, 1.71385, 0.299598, 0.000505937, 0.0647114, 5.01559e-6, 1.49885e-5, 0.00544226, 4.35115e-6, 3.03577e-5, 2.15171e-5, 13.9863, 2.77415, 3.06754, 3.09629, 0.129791, 0.002272, 2.12184, 0.461891, 0.362844, 0.0293403, 0.000990191, 0.00063532, 2.98906e-5, 3.40495e-6, 0.000103298, 5.23445e-7, 7.58079e-6, 0.00361932, 1.11528, 1.93197, 3.82999, 3.8369, 7.05901e-6, 3.52973e-5, 0.714193, 0.0650093, 0.521932, 0.0465728, 0.00111718, 6.65118e-6, 0.00254141, 0.0378996, 9.51571e-5, 7.83845e-7, 0.000242097, 0.00102866, 1.5658, 3.63417, 3.88274, 3.97669, 0.0956175, 0.64807, 2.13597, 0.396049, 1.0315, 0.395203, 0.00162117, 0.00053917, 0.000125216, 0.0929912, 0.00102563, 0.000176756, 0.00184695, 0.0261891, 13.4388, 0.395591, 0.0829676, 0.000898367, 0.0907073, 0.000483117, 0.0306672, 0.166496, 0.000362649, 0.399909, 0.00273783, 0.0227258, 3.30817e-6, 8.66783e-6, 9.05333e-6, 2.4177e-5, 0.000184217, 2.80325e-6, 3.56209, 1.07825, 3.95194, 0.447834, 0.0389113, 0.110297, 2.13567, 0.27307, 3.36894, 0.368404, 0.000296116, 0.0415469, 4.95713e-6, 2.47911e-5, 2.06936e-5, 3.39401e-5, 0.0908003, 3.87974e-6, 9.97915, 9.99163, 0.288515, 0.632015, 2.12838, 0.320644, 6.53941, 0.325873]

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

    x = ["BEZ235", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib", "Trametinib"]
    dfs1 = DataFrame(x=x, y=gi[1, 1, :], y2=gi[1, 2, :], label="G1")
    dfs2 = DataFrame(x=x, y=gi[2, 1, :], y2=gi[2, 2, :], label="SG2")
    dfs3 = DataFrame(x=x, y=sum(deathContG1, dims = 1)[1, :], y2=sum(deathEC50G1, dims = 1)[1, :], label="G1")
    dfs4 = DataFrame(x=x, y=sum(deathContG2, dims = 1)[1, :], y2=sum(deathEC50G2, dims = 1)[1, :], label="SG2")

    p = [Gadfly.plot(dfs1, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=0.0,ymax=60.0), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]"), Coord.Cartesian(ymin=0.0,ymax=60.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"])),
    
    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=-0.05,ymax=0.2), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"), Coord.Cartesian(ymin=-0.05,ymax=0.2), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"]))]
    
    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("figure81.svg", 8inch, 8inch), pl1)
end