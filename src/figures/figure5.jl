""" Plot the 6 drugs from HCC dataset. """
function mt21_single(tensr, condition, g1g2, titles=false)
    t = LinRange(0.0, 95, size(tensr)[1])
    dfs = [DataFrame(x=t, y=tensr[:, i], conc=condition[i]) for i = 1:8]
    DF = vcat(dfs...)
    if titles == false
        return Gadfly.plot(DF, x="x", y="y", color="conc", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"))
    else
        return Gadfly.plot(DF, x="x", y="y", color="conc", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"), Guide.title(titles))        
    end
end

function figure5()
    setGadflyTheme()

    tens, names, _, conds, _, _ = DrugResponseModel.hcc_all()

    pp = []
    for i = 1:size(tens)[4]
        push!(pp, DrugResponseModel.Eachdrug(tens[1, :, :, i], conds[i], "G1", names[i]))
        push!(pp, DrugResponseModel.Eachdrug(tens[2, :, :, i], conds[i], "S/G2", names[i]))
    end

    pl = plotGrid((3, 4), [pp...];)
    return draw(SVG("figure5.svg", 16inch, 12inch), pl)
end

function figure51()
    setGadflyTheme()

    tens, names, _, conds, _, _ = DrugResponseModel.mt1_all()

    pp = []
    for i = 1:size(tens)[4]
        push!(pp, DrugResponseModel.mt21_single(tens[1, :, :, i], conds[i], "G1", names[i]))
        push!(pp, DrugResponseModel.mt21_single(tens[2, :, :, i], conds[i], "S/G2", names[i]))
    end

    pl = plotGrid((3, 4), [pp...];)
    return draw(SVG("figure52.svg", 16inch, 12inch), pl)
end

function figure54()
    setGadflyTheme()

    tens, names, _, conds, _, _ = DrugResponseModel.mt1_all()

    pp = []
    for i = 1:size(tens)[4]
        push!(pp, DrugResponseModel.mt21_single(tens[1, :, :, i], conds[i], "G1", names[i]))
        push!(pp, DrugResponseModel.mt21_single(tens[2, :, :, i], conds[i], "S/G2", names[i]))
    end

    pl = plotGrid((3, 4), [pp...];)
    return draw(SVG("figure53.svg", 16inch, 12inch), pl)
end