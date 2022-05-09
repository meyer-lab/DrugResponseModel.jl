""" Figure 8: the quantified effects of drugs for HCC cell line when fitting at once. """
function plt_2(i, j, eff_name)
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end

    params = [4.83332, 0.613402, 2.74895e-5, 0.00959796, 0.898324, 0.319393, 0.431104, 1.8271, 1.80044, 0.00417866, 0.0477392, 1.86715e-6, 3.26458e-6, 1.91951e-6, 4.58963e-7, 1.86302e-6, 4.62873e-7, 1.74927e-7, 173.145, 9.9888, 2.42259e-5, 0.304448, 0.433302, 0.0559058, 0.0241491, 1.80361, 1.77422, 0.0563967, 7.53949e-6, 3.69925, 0.000105467, 1.49994e-6, 5.84926e-7, 7.08798e-5, 1.34358e-5, 8.45627e-7, 0.857942, 2.65903, 0.158182, 0.0277238, 3.01433e-6, 0.0329552, 0.0170167, 1.81469, 1.78454, 0.00761463, 2.22794e-5, 1.05081e-6, 6.02528e-6, 3.91612e-6, 0.000520593, 2.96725e-5, 2.11235e-5, 0.00010379, 1.02373, 9.99138, 1.27734, 0.298102, 1.67355, 0.446294, 0.180929, 1.81345, 1.78349, 0.325204, 6.71652e-6, 4.49531e-6, 0.404435, 1.01417e-7, 6.46561e-7, 7.86395e-6, 2.4704e-6, 1.8778e-7, 11.9208, 0.245903, 0.233606, 0.0729809, 0.58465, 0.395487, 0.280238, 1.81472, 1.7865, 6.30423e-6, 1.40879e-6, 1.07929e-6, 2.12691e-6, 5.48223e-7, 2.71e-7, 1.24475e-7, 1.61646e-6, 7.1574e-8, 21.3821, 1.12372, 0.249267, 0.153966, 1.42807, 0.328861, 0.120697, 1.80585, 0.474918, 0.0569997, 2.13323e-5, 1.02955e-6, 9.57608e-6, 5.97748e-7, 1.35781e-7, 5.21108e-6, 7.20512e-7, 0.000360178, 1.27936, 0.295626, 1.62929, 0.114501, 0.193466, 1.81144, 1.78166, 0.293]
    pp = getODEparams(params, cs)

    dfs = [DataFrame(x=concs[i], y=pp[Int(j+k-1), :, i], label=string(eff_name*"$k")) for k = 1:4]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.title(names[i]), Guide.xlabel("concentrations [nM]"), Guide.ylabel(eff_name), style(line_width=1mm), Coord.Cartesian(ymin=0.0,ymax=4.0))
end

function figure8()
    setGadflyTheme()
    plts = []
    for i=1:6
        push!(plts, plt_2(i, 1, "G1 progression"))
        push!(plts, plt_2(i, 5, "S-G2 progression"))
        push!(plts, plt_2(i, 9, "G1 cell death"))
        push!(plts, plt_2(i, 12, "S-G2 cell death"))
    end

    pl1 = plotGrid((6, 4), [plts...];)
    return draw(SVG("figure8.svg", 16inch, 20inch), pl1)
end