""" Figure 8: the quantified effects of drugs for HCC cell line when fitting at once. """

""" Helper function for plotting the effects. """
function effectsplot_helper(params, cs)
    efcs = getODEparams(params, cs)
    # params at EC50
    ec50 = zeros(16, 4)
    conc_ec50 = zeros((1, 4))
    conc_ec50[1, 1:4] = [params[18*k+1] for k = 0:3]
    ec50 = getODEparams(params, conc_ec50)[:, 1, :]
    print(conc_ec50)

    # phase durations
    # @ control
    gi = zeros(2, 2, 4)
    gi[1, 1, :] .= (2 ./ efcs[1, 1, :] .+ 2 ./ efcs[2, 1, :] .+ 2 ./ efcs[3, 1, :] .+ 2 ./ efcs[4, 1, :])
    gi[2, 1, :] .= (5 ./ efcs[5, 1, :] .+ 5 ./ efcs[6, 1, :] .+ 5 ./ efcs[7, 1, :] .+ 5 ./ efcs[8, 1, :])

    # @ ec50
    gi[1, 2, :] .= (2 ./ ec50[1, :] .+ 2 ./ ec50[2, :] .+ 2 ./ ec50[3, :] .+ 2 ./ ec50[4, :])
    gi[2, 2, :] .= (5 ./ ec50[5, :] .+ 5 ./ ec50[6, :] .+ 5 ./ ec50[7, :] .+ 5 ./ ec50[8, :])

    # cell deaths
    deathContG1 = zeros(4, 4)
    deathEC50G1 = zeros(4, 4)
    deathContG1[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[1, 1, :])
    deathEC50G1[1, :] .= (ec50[9, :]) ./ (ec50[9, :] .+ ec50[1, :])
    for i = 2:4
        deathContG1[i, :] .= (1 .- deathContG1[i - 1, :]) .* (efcs[i + 8, 1, :]) ./ (efcs[i + 8, 1, :] .+ efcs[i, 1, :])
        deathEC50G1[i, :] .= (1 .- deathEC50G1[i - 1, :]) .* (ec50[i + 8, :]) ./ (ec50[i + 8, :] .+ ec50[i, :])
    end
    deathContG2 = zeros(4, 4)
    deathEC50G2 = zeros(4, 4)
    deathContG2[1, :] .= (efcs[13, 1, :]) ./ (efcs[13, 1, :] .+ efcs[5, 1, :])
    deathEC50G2[1, :] .= (ec50[13, :]) ./ (ec50[13, :] .+ ec50[5, :])
    for i = 14:16
        deathContG2[i - 12, :] = (1 .- deathContG2[i - 13, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 8, 1, :])
        deathEC50G2[i - 12, :] = (1 .- deathEC50G2[i - 13, :]) .* (ec50[i, :]) ./ (ec50[i, :] .+ ec50[i - 8, :])
    end

    return gi, deathContG1, deathEC50G1, deathContG2, deathEC50G2

end

""" HCC1143 effect plot """
function figure81()
    tensor, names, concs, conds = DrugResponseModel.hcc_all()
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end

    params = [1.36845, 4.05584, 3.99669, 3.99714, 0.102503, 0.458295, 0.659276, 3.99775, 0.365036, 0.627148, 6.94891e-5, 0.000136151, 0.0127684, 0.0517801, 3.20979e-5, 4.19056e-5, 7.37436e-6, 0.0228549, 1.11308, 0.264787, 3.99759, 1.45882e-5, 3.06579e-6, 0.375592, 0.102024, 3.99633, 0.320449, 0.722804, 3.53551e-6, 0.032066, 1.67137e-7, 4.27832e-7, 8.83771e-8, 2.89182e-6, 4.5214e-7, 6.8609e-7, 16.4934, 2.45814, 1.12299, 1.55116, 0.0033046, 0.185649, 0.673943, 0.305188, 0.0594824, 0.0924852, 5.03661e-5, 7.25378e-6, 0.00764243, 3.2136e-6, 1.2219e-6, 7.56255e-7, 0.00440004, 0.00037535, 0.741965, 1.89296, 3.98903, 3.99547, 2.31225e-6, 0.383961, 0.464131, 2.02298, 0.091986, 0.293161, 1.52185e-5, 2.38363e-5, 0.021285, 1.42916e-6, 3.24183e-8, 3.92282e-6, 1.48194e-7, 1.03137e-6, 3.99914, 3.99967, 0.482352, 0.358652, 0.531748, 3.99842, 0.35508, 0.537322]

    gi, deathContG1, deathEC50G1, deathContG2, deathEC50G2 = effectsplot_helper(params, cs)

    x = ["Paclitaxel", "Palbociclib", "Doxorubicin", "Gemcitabine"]
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


""" 21MT1 effect plot """
function figure82()
    ENV["GKSwstype"]="nul"
    ten, nams, concs, conds = DrugResponseModel.mt1_all()

    cs = zeros(8, 6)
    for i=1:6
        cs[:, i] .= concs[i]
    end

    params = [1.57556, 1.62607, 4.35979e-6, 2.20023, 6.40848e-9, 1.63139e-7, 4.0, 1.03588, 2.82295e-7, 1.518, 7.26753e-9, 3.32831e-8, 7.38826e-9, 0.0847927, 1.59243e-8, 6.84661e-9, 0.0635307, 2.0618e-9, 18.4637, 0.498753, 3.99999, 4.0, 0.00330672, 0.509792, 4.0, 0.669307, 3.92732, 1.30048, 9.84433e-8, 1.17613e-8, 1.03609e-9, 1.63366e-8, 4.43533e-8, 1.11375e-9, 1.15564e-8, 2.01859e-9, 41.7844, 1.70309, 2.01363e-6, 0.0480213, 1.40979e-9, 1.5706e-7, 3.99999, 0.89413, 4.0, 7.42191e-9, 1.00757e-6, 0.0640029, 7.28897e-9, 1.06527, 2.77222e-8, 4.83221e-9, 4.92092e-9, 1.18316e-9, 0.733685, 1.82624, 4.0, 0.0108699, 0.991755, 0.0497405, 3.99663, 0.614765, 2.77109e-9, 1.10172, 1.66923e-7, 0.00238703, 6.67847e-8, 1.01807e-9, 5.22232e-9, 2.70425e-9, 1.80354e-9, 2.04378e-9, 4.0, 4.0, 0.574277, 2.19519, 3.9953, 0.627565, 4.0, 1.23897]

    gi, deathContG1, deathEC50G1, deathContG2, deathEC50G2 = effectsplot_helper(params, cs)

    x = ["Paclitaxel", "Palbociclib", "Doxorubicin", "Gemcitabine"]
    dfs1 = DataFrame(x=x, y=gi[1, 1, :], y2=gi[1, 2, :], label="G1")
    dfs2 = DataFrame(x=x, y=gi[2, 1, :], y2=gi[2, 2, :], label="SG2")
    dfs3 = DataFrame(x=x, y=sum(deathContG1, dims = 1)[1, :], y2=sum(deathEC50G1, dims = 1)[1, :], label="G1")
    dfs4 = DataFrame(x=x, y=sum(deathContG2, dims = 1)[1, :], y2=sum(deathEC50G2, dims = 1)[1, :], label="SG2")

    p = [Gadfly.plot(dfs1, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=0.0,ymax=25.0), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]"), Coord.Cartesian(ymin=0.0,ymax=25.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"])),
    
    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=-0.05,ymax=1.0), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"), Coord.Cartesian(ymin=-0.05,ymax=1.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"]))]

    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("figure82.svg", 8inch, 8inch), pl1)
end


""" MDAMB-157 effect plot """
function figure83()
    ENV["GKSwstype"]="nul"
    ten, nams, concs, conds = DrugResponseModel.mda_all()

    indexes = [1, 2, 7, 8, 9, 10]
    cs = zeros(8, 6)
    for i=1:6
        cs[:, i] .= concs[indexes[i]]
    end

    params = [2.13839, 2.83884, 1.51441, 1.50356, 0.0863139, 0.289618, 0.753405, 7.01123e-6, 0.0869785, 3.99967, 1.53553e-5, 7.3733e-6, 1.41833e-7, 0.0528917, 0.0221832, 0.374224, 2.10944e-7, 7.27666e-8, 412.122, 0.685044, 5.93963e-6, 1.44106e-6, 0.269538, 0.0777158, 0.345028, 3.99898, 0.305618, 3.99972, 0.0880853, 0.0556922, 5.90658e-7, 6.21931e-7, 6.23376e-8, 3.69143e-7, 3.19237e-8, 4.25434e-6, 12.8393, 1.10132, 4.33747e-7, 3.75999e-7, 0.258907, 4.06523e-6, 3.42362e-7, 0.190433, 0.267539, 0.18205, 0.0253506, 0.0219874, 0.00864397, 0.00222323, 0.00597668, 4.48199e-7, 1.28088e-6, 1.45763e-7, 0.471157, 5.8535, 2.15688, 1.73648, 0.253272, 0.13052, 0.160907, 0.129681, 0.464879, 3.99994, 0.0109992, 0.0343582, 9.64949e-7, 9.13167e-8, 1.38366e-8, 0.01639, 0.0109108, 0.000553012, 3.99994, 3.99988, 0.104727, 1.34732, 0.475374, 3.99997, 0.298439, 3.99996]

    gi, deathContG1, deathEC50G1, deathContG2, deathEC50G2 = effectsplot_helper(params, cs)

    x = ["Paclitaxel", "Palbociclib", "Doxorubicin", "Gemcitabine"]
    dfs1 = DataFrame(x=x, y=gi[1, 1, :], y2=gi[1, 2, :], label="G1")
    dfs2 = DataFrame(x=x, y=gi[2, 1, :], y2=gi[2, 2, :], label="SG2")
    dfs3 = DataFrame(x=x, y=sum(deathContG1, dims = 1)[1, :], y2=sum(deathEC50G1, dims = 1)[1, :], label="G1")
    dfs4 = DataFrame(x=x, y=sum(deathContG2, dims = 1)[1, :], y2=sum(deathEC50G2, dims = 1)[1, :], label="SG2")

    p = [Gadfly.plot(dfs1, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=0.0,ymax=50.0), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]"), Coord.Cartesian(ymin=0.0,ymax=50.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"])),
    
    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.manual_color_key("", ["G1_control", "G1_ec50"], ["black", "red"]), Coord.Cartesian(ymin=-0.05,ymax=1.0), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),
    
    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, Theme(default_color=colorant"black")), 
    layer(x="x", y="y2", Geom.point, Theme(default_color=colorant"red")),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"), Coord.Cartesian(ymin=-0.05,ymax=1.0), Guide.manual_color_key("", ["SG2_control", "SG2_ec50"], ["black", "red"]))]

    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("figure83.svg", 8inch, 8inch), pl1)
end

""" Plotting the effects of all three cell lines on top of each other on a plot. """
function figure84()

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

    params_hcc = [1.36845, 4.05584, 3.99669, 3.99714, 0.102503, 0.458295, 0.659276, 3.99775, 0.365036, 0.627148, 6.94891e-5, 0.000136151, 0.0127684, 0.0517801, 3.20979e-5, 4.19056e-5, 7.37436e-6, 0.0228549, 1.11308, 0.264787, 3.99759, 1.45882e-5, 3.06579e-6, 0.375592, 0.102024, 3.99633, 0.320449, 0.722804, 3.53551e-6, 0.032066, 1.67137e-7, 4.27832e-7, 8.83771e-8, 2.89182e-6, 4.5214e-7, 6.8609e-7, 16.4934, 2.45814, 1.12299, 1.55116, 0.0033046, 0.185649, 0.673943, 0.305188, 0.0594824, 0.0924852, 5.03661e-5, 7.25378e-6, 0.00764243, 3.2136e-6, 1.2219e-6, 7.56255e-7, 0.00440004, 0.00037535, 0.741965, 1.89296, 3.98903, 3.99547, 2.31225e-6, 0.383961, 0.464131, 2.02298, 0.091986, 0.293161, 1.52185e-5, 2.38363e-5, 0.021285, 1.42916e-6, 3.24183e-8, 3.92282e-6, 1.48194e-7, 1.03137e-6, 3.99914, 3.99967, 0.482352, 0.358652, 0.531748, 3.99842, 0.35508, 0.537322]
    params_21mt1 = [1.57556, 1.62607, 4.35979e-6, 2.20023, 6.40848e-9, 1.63139e-7, 4.0, 1.03588, 2.82295e-7, 1.518, 7.26753e-9, 3.32831e-8, 7.38826e-9, 0.0847927, 1.59243e-8, 6.84661e-9, 0.0635307, 2.0618e-9, 18.4637, 0.498753, 3.99999, 4.0, 0.00330672, 0.509792, 4.0, 0.669307, 3.92732, 1.30048, 9.84433e-8, 1.17613e-8, 1.03609e-9, 1.63366e-8, 4.43533e-8, 1.11375e-9, 1.15564e-8, 2.01859e-9, 41.7844, 1.70309, 2.01363e-6, 0.0480213, 1.40979e-9, 1.5706e-7, 3.99999, 0.89413, 4.0, 7.42191e-9, 1.00757e-6, 0.0640029, 7.28897e-9, 1.06527, 2.77222e-8, 4.83221e-9, 4.92092e-9, 1.18316e-9, 0.733685, 1.82624, 4.0, 0.0108699, 0.991755, 0.0497405, 3.99663, 0.614765, 2.77109e-9, 1.10172, 1.66923e-7, 0.00238703, 6.67847e-8, 1.01807e-9, 5.22232e-9, 2.70425e-9, 1.80354e-9, 2.04378e-9, 4.0, 4.0, 0.574277, 2.19519, 3.9953, 0.627565, 4.0, 1.23897]
    params_mdamb = [2.13839, 2.83884, 1.51441, 1.50356, 0.0863139, 0.289618, 0.753405, 7.01123e-6, 0.0869785, 3.99967, 1.53553e-5, 7.3733e-6, 1.41833e-7, 0.0528917, 0.0221832, 0.374224, 2.10944e-7, 7.27666e-8, 412.122, 0.685044, 5.93963e-6, 1.44106e-6, 0.269538, 0.0777158, 0.345028, 3.99898, 0.305618, 3.99972, 0.0880853, 0.0556922, 5.90658e-7, 6.21931e-7, 6.23376e-8, 3.69143e-7, 3.19237e-8, 4.25434e-6, 12.8393, 1.10132, 4.33747e-7, 3.75999e-7, 0.258907, 4.06523e-6, 3.42362e-7, 0.190433, 0.267539, 0.18205, 0.0253506, 0.0219874, 0.00864397, 0.00222323, 0.00597668, 4.48199e-7, 1.28088e-6, 1.45763e-7, 0.471157, 5.8535, 2.15688, 1.73648, 0.253272, 0.13052, 0.160907, 0.129681, 0.464879, 3.99994, 0.0109992, 0.0343582, 9.64949e-7, 9.13167e-8, 1.38366e-8, 0.01639, 0.0109108, 0.000553012, 3.99994, 3.99988, 0.104727, 1.34732, 0.475374, 3.99997, 0.298439, 3.99996]
    
    gi_h, deathContG1_h, deathEC50G1_h, deathContG2_h, deathEC50G2_h = effectsplot_helper(params_hcc, cs[:, :, 1])
    gi_t, deathContG1_t, deathEC50G1_t, deathContG2_t, deathEC50G2_t = effectsplot_helper(params_21mt1, cs[:, :, 2])
    gi_m, deathContG1_m, deathEC50G1_m, deathContG2_m, deathEC50G2_m = effectsplot_helper(params_mdamb, cs[:, :, 3])

    x = ["Paclitaxel", "Palbociclib", "Doxorubicin", "Gemcitabine"]
    dfs1 = DataFrame(x=repeat(x, 3), y=vcat(gi_h[1, 1, :], gi_t[1, 1, :], gi_m[1, 1, :]), y2=vcat(gi_h[1, 2, :], gi_t[1, 2, :], gi_m[1, 2, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=4), label="G1")
    dfs2 = DataFrame(x=repeat(x, 3), y=vcat(gi_h[2, 1, :], gi_t[2, 1, :], gi_m[2, 1, :]), y2=vcat(gi_h[2, 2, :], gi_t[2, 2, :], gi_m[2, 2, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=4), label="SG2")

    dfs3 = DataFrame(x=repeat(x, 3), y=vcat(sum(deathContG1_h, dims = 1)[1, :], sum(deathContG1_t, dims = 1)[1, :], sum(deathContG1_m, dims = 1)[1, :]), y2=vcat(sum(deathEC50G1_h, dims = 1)[1, :], sum(deathEC50G1_t, dims = 1)[1, :], sum(deathEC50G1_m, dims = 1)[1, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=4), label="G1")
    dfs4 = DataFrame(x=repeat(x, 3), y=vcat(sum(deathContG2_h, dims = 1)[1, :], sum(deathContG2_t, dims = 1)[1, :], sum(deathContG2_m, dims = 1)[1, :]), y2=vcat(sum(deathEC50G2_h, dims = 1)[1, :], sum(deathEC50G2_t, dims = 1)[1, :], sum(deathEC50G2_m, dims = 1)[1, :]), cellline=repeat(["HCC1143", "21mt1", "MDAMB157"], inner=4), label="SG2")

    p = [Gadfly.plot(dfs1, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Coord.Cartesian(ymin=0.0,ymax=50.0), Guide.ylabel("phase duration [hr]"), Guide.title("G1 drug effects")),

    Gadfly.plot(dfs2, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("phase duration [hr]"), Coord.Cartesian(ymin=0.0,ymax=50.0)),

    Gadfly.plot(dfs3, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.ylabel("cell death probability"), Guide.title("G1 drug effects")),

    Gadfly.plot(dfs4, layer(x="x", y="y", Geom.point, size=[0.07sx], shape=[Gadfly.Shape.cross], color="cellline"), 
    layer(x="x", y="y2", Geom.point, shape=[Gadfly.Shape.hexagon], size=[0.07sx], color="cellline"),
    Guide.xlabel("drugs"), Guide.title("SG2 drug effects"), Guide.ylabel("cell death probability"))]

    pl1 = plotGrid((2, 2), [p...];)
    return draw(SVG("figure84.svg", 10inch, 6inch), pl1)
end