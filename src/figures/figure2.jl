""" Plotting the model and data. """

drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

function Eachdrug_sim(G_sim, G_data, condition, g1g2, drug_name)
    t = LinRange(0.0, 95, 189)
    df = [DataFrame(x=t, y=G_sim[:, i], y2=G_data[:, i], conc=condition[i]) for i = 1:8]
    DF = vcat(df...)
    Gadfly.plot(DF, layer(x="x", y="y", color="conc", Geom.line), 
                    layer(x="x", y="y2", color="conc", Geom.line, Theme(line_style=[:dash])), 
                    Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"), Coord.Cartesian(ymin=-0.05,ymax=2.0), Guide.title(drug_name))
end

function figure2()
    setGadflyTheme()

    t = LinRange(0.0, 95.5, 189)
    plots = []
    for drug in drugs
        g1, g2, conc = DrugResponseModel.load_newData(drug)
        fitness, p = DrugResponseModel.optimize_hill(conc, g1[:, :, 1], g2[:, :, 1])
        pODE = DrugResponseModel.getODEparams(p, conc)
        Gs = zeros(189, 8, 2)
        for i=1:8
            Gs[:, i, 1], Gs[:, i, 2], _ = DrugResponseModel.predict(pODE[:, i, 1], pODE[:, 1, 1], t)
        end
        push!(plots, Eachdrug_sim(Gs[:, :, 1], g1[:, :, 1], conc, "G1", drug), Eachdrug_sim(Gs[:, :, 2], g2[:, :, 1], conc, "S/G2", drug))
    end

    pl = plotGrid((4, 6), [plots..., nothing, nothing])
    return draw(SVG("figure2.svg", 22inch, 18inch), pl)
end