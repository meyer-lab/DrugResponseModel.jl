""" plot the new data. """

function Eachdrug(tensr, condition, g1g2)
    t = LinRange(0.0, 95, 189)
    dfs = [DataFrame(x=t, y=tensr[:, i], conc=condition[i]) for i = 1:8]
    DF = vcat(dfs...)
    return Gadfly.plot(DF, x="x", y="y", color="conc", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"), Coord.Cartesian(ymin=-0.05,ymax=2.0))
end

function figure1()
    setGadflyTheme()
    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)[:, 1, :, :, :]

    pp = []
    for i = 1:11
        push!(pp, Eachdrug(tens[1, :, :, i], conditions[i], "G1"))
        push!(pp, Eachdrug(tens[2, :, :, i], conditions[i], "S/G2"))
    end

    pl = plotGrid((6, 4), [pp..., nothing, nothing];)
    return draw(SVG("figure1.svg", 18inch, 20inch), pl)
end

function figure2()
    setGadflyTheme()
    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)
    pp = []
    for i = 1:11
        push!(pp, Eachdrug(tens[1, 1, :, :, i] .+ tens[2, 1, :, :, i], conditions[i], "total"))
    end

    pl = plotGrid((3, 4), [pp..., nothing];)
    return draw(SVG("figure2.svg", 18inch, 14inch), pl)
end
