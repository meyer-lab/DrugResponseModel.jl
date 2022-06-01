""" Figure 8: the quantified effects of drugs for HCC cell line when fitting at once. """
function plt_2(i, j, eff_name, ymax)
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end

    params = [6.76402, 1.45485, 0.0901376, 2.00322, 2.00283, 0.187222, 2.005, 2.00412, 0.148397, 0.0586825, 0.0213053, 2.52195e-7, 2.41545e-8, 1.44707e-8, 3.29222e-8, 1.37208e-8, 5.59402e-9, 0.00126715, 43.3336, 1.39287, 0.139558, 0.138035, 0.0904822, 0.0493381, 9.62971e-8, 5.98781e-7, 0.264581, 0.00759821, 6.58923e-9, 2.90942e-8, 7.08595e-8, 4.87951e-8, 2.21321e-8, 1.84993e-8, 0.0515671, 6.56555e-9, 1.05139, 2.02354, 0.0172485, 2.00418, 0.258644, 0.949032, 2.002, 2.00531, 1.1072e-7, 0.119229, 2.07555e-8, 0.407274, 8.36152e-8, 9.54029e-8, 2.32011e-7, 4.88497e-8, 0.0021182, 4.3538e-8, 1.14876, 3.86715, 0.15917, 2.01371, 2.01491, 1.44785, 2.00591, 2.01811, 0.140134, 0.749386, 0.0433281, 1.6227e-7, 5.84897e-8, 5.88594e-8, 2.78218e-8, 1.63851e-8, 5.49249e-9, 6.74174e-8, 5000.0, 0.310183, 0.172035, 2.00422, 2.00337, 1.4227, 2.00219, 2.00051, 0.0141024, 0.236535, 0.0479571, 2.37249e-7, 1.46444e-7, 5.86096e-8, 2.4064e-8, 1.12806e-8, 7.34672e-9, 8.34734e-8, 35.3323, 1.03572, 0.162309, 2.00718, 2.00352, 1.41681, 2.00072, 2.0004, 0.0604097, 0.662226, 1.67197e-8, 1.02143e-7, 5.22464e-8, 2.72402e-8, 2.40392e-8, 7.06357e-9, 8.59539e-9, 0.0646765, 0.0750058, 2.0, 1.99998, 1.41318, 2.0, 2.0, 0.126361, 0.810005]
    pp = getODEparams(params, cs)

    dfs = [DataFrame(x=concs[i], y=pp[Int(j+k-1), :, i], label=string(eff_name*"$k")) for k = 1:4]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.title(names[i]), Guide.xlabel("concentrations [nM]"), Guide.ylabel(eff_name), style(line_width=1mm), Coord.Cartesian(ymin=0.0,ymax=ymax))
end

function figure8()
    setGadflyTheme()
    plts = []
    for i=1:6
        push!(plts, plt_2(i, 1, "G1 progression", 2.0))
        push!(plts, plt_2(i, 5, "S-G2 progression", 2.0))
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

    params = [6.76402, 1.45485, 0.0901376, 2.00322, 2.00283, 0.187222, 2.005, 2.00412, 0.148397, 0.0586825, 0.0213053, 2.52195e-7, 2.41545e-8, 1.44707e-8, 3.29222e-8, 1.37208e-8, 5.59402e-9, 0.00126715, 43.3336, 1.39287, 0.139558, 0.138035, 0.0904822, 0.0493381, 9.62971e-8, 5.98781e-7, 0.264581, 0.00759821, 6.58923e-9, 2.90942e-8, 7.08595e-8, 4.87951e-8, 2.21321e-8, 1.84993e-8, 0.0515671, 6.56555e-9, 1.05139, 2.02354, 0.0172485, 2.00418, 0.258644, 0.949032, 2.002, 2.00531, 1.1072e-7, 0.119229, 2.07555e-8, 0.407274, 8.36152e-8, 9.54029e-8, 2.32011e-7, 4.88497e-8, 0.0021182, 4.3538e-8, 1.14876, 3.86715, 0.15917, 2.01371, 2.01491, 1.44785, 2.00591, 2.01811, 0.140134, 0.749386, 0.0433281, 1.6227e-7, 5.84897e-8, 5.88594e-8, 2.78218e-8, 1.63851e-8, 5.49249e-9, 6.74174e-8, 5000.0, 0.310183, 0.172035, 2.00422, 2.00337, 1.4227, 2.00219, 2.00051, 0.0141024, 0.236535, 0.0479571, 2.37249e-7, 1.46444e-7, 5.86096e-8, 2.4064e-8, 1.12806e-8, 7.34672e-9, 8.34734e-8, 35.3323, 1.03572, 0.162309, 2.00718, 2.00352, 1.41681, 2.00072, 2.0004, 0.0604097, 0.662226, 1.67197e-8, 1.02143e-7, 5.22464e-8, 2.72402e-8, 2.40392e-8, 7.06357e-9, 8.59539e-9, 0.0646765, 0.0750058, 2.0, 1.99998, 1.41318, 2.0, 2.0, 0.126361, 0.810005]

    efcs = getODEparams(params, cs)
    # params at EC50
    ec50 = zeros(16, 6)
    conc_ec50 = zeros((1, 6))
    conc_ec50[1, 1:6] = [params[18*k+1] for k = 0:5]
    ec50 = getODEparams(params, conc_ec50)[:, 1, :]

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
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=0.0,ymax=90.0), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]"), Coord.Cartesian(ymin=0.0,ymax=90.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"])),
    
    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=-0.05,ymax=0.25), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"), Coord.Cartesian(ymin=-0.05,ymax=0.25), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"]))]
    
    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("figure81.svg", 8inch, 8inch), pl1)
end