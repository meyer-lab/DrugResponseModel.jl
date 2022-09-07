""" Figure4 plots the effects over concentrations. """
drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

function plt_1(i, j, eff_name)
    pp = DrugResponseModel.out_ODEparams() # 11 drugs x 16 parameters x 8 concentrations
    _, _, conc = DrugResponseModel.load_newData(drugs[i])
    cs = convert(Array{Float64,1}, conc)
    dfs = [DataFrame(x=cs, y=pp[i, Int(j+k-1), :], label=string(eff_name*"$k")) for k = 1:4]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="label", Geom.line, Guide.title(drugs[i]), Guide.xlabel("concentrations [nM]"), Guide.ylabel(eff_name), Coord.Cartesian(ymin=0.0,ymax=10.0))
end

function figure4()
    setGadflyTheme()
    plts = []
    for i=1:6
        push!(plts, plt_1(i, 1, "α"))
        push!(plts, plt_1(i, 5, "β"))
        push!(plts, plt_1(i, 9, "γ_1"))
        push!(plts, plt_1(i, 12, "γ_2"))
    end

    pl1 = plotGrid((6, 4), [plts...];)
    return draw(SVG("figure4a.svg", 12inch, 20inch), pl1)
end

function figure41()
    setGadflyTheme()
    plts = []
    for i=7:11
        push!(plts, plt_1(i, 1, "α"))
        push!(plts, plt_1(i, 5, "β"))
        push!(plts, plt_1(i, 9, "γ_1"))
        push!(plts, plt_1(i, 12, "γ_2"))
    end

    pl1 = DrugResponseModel.plotGrid((5, 4), [plts...];)
    return draw(SVG("figure4b.svg", 12inch, 18inch), pl1)
end
