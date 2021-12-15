""" Figure4 plots the effects over concentrations. """
drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

function plt_1(i, j, eff_name)
    pp = DrugResponseModel.out_ODEparams() # 11 drugs x 16 parameters x 8 concentrations
    _, _, conc = DrugResponseModel.load_newData(drugs[i])
    dfs = [DataFrame(x=conc, y=pp[i, j+k-1, :], label=string(eff_name*"$k")) for k = 1:4]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.xlabel("concentrations [nM]"), Guide.ylabel(eff_name), Coord.Cartesian(ymin=0.0,ymax=10.0))
end

function figure4()
    setGadflyTheme()
    plts = []
    for (i, drug) in enumerate(drugs[1:8])
        push!(plts, plt_1(i, 1, "α"))
        push!(plts, plt_1(i, 5, "β"))
        push!(plts, plt_1(i, 9, "γ_1"))
        push!(plts, plt_1(i, 12, "γ_2"))
    end

    pl1 = DrugResponseModel.plotGrid((8, 4), [plts...];)
    return draw(SVG("figure4.svg", 4inch, 14inch), pl1)
end

function figure5()
    setGadflyTheme()
    plts = []
    for (i, drug) in enumerate(drugs[9:11])
        push!(plts, plt_1(i, 1, "α"))
        push!(plts, plt_1(i, 5, "β"))
        push!(plts, plt_1(i, 9, "γ_1"))
        push!(plts, plt_1(i, 12, "γ_2"))
    end

    pl1 = DrugResponseModel.plotGrid((3, 4), [plts...];)
    return draw(SVG("figure5.svg", 10inch, 14inch), pl1)
end
