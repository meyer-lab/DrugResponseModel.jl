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
    tens = mean(tensor, dims=2)

    p1 = Eachdrug(tens[1, 1, :, :, 1], conditions[1], "G1")
    p2 = Eachdrug(tens[2, 1, :, :, 1], conditions[1], "S/G2")
    
    p3 = Eachdrug(tens[1, 1, :, :, 2], conditions[2], "G1")
    p4 = Eachdrug(tens[2, 1, :, :, 2], conditions[2], "S/G2")

    p5 = Eachdrug(tens[1, 1, :, :, 3], conditions[3], "G1")
    p6 = Eachdrug(tens[2, 1, :, :, 3], conditions[3], "S/G2")

    p7 = Eachdrug(tens[1, 1, :, :, 4], conditions[4], "G1")
    p8 = Eachdrug(tens[2, 1, :, :, 4], conditions[4], "S/G2")

    p9 = Eachdrug(tens[1, 1, :, :, 5], conditions[5], "G1")
    p10 = Eachdrug(tens[2, 1, :, :, 5], conditions[5], "S/G2")

    p11 = Eachdrug(tens[1, 1, :, :, 6], conditions[6], "G1")
    p12 = Eachdrug(tens[2, 1, :, :, 6], conditions[6], "S/G2")

    p13 = Eachdrug(tens[1, 1, :, :, 7], conditions[7], "G1")
    p14 = Eachdrug(tens[2, 1, :, :, 7], conditions[7], "S/G2")
    
    p15 = Eachdrug(tens[1, 1, :, :, 8], conditions[8], "G1")
    p16 = Eachdrug(tens[2, 1, :, :, 8], conditions[8], "S/G2")

    p17 = Eachdrug(tens[1, 1, :, :, 9], conditions[9], "G1")
    p18 = Eachdrug(tens[2, 1, :, :, 9], conditions[9], "S/G2")

    p19 = Eachdrug(tens[1, 1, :, :, 10], conditions[10], "G1")
    p20 = Eachdrug(tens[2, 1, :, :, 10], conditions[10], "S/G2")

    p21 = Eachdrug(tens[1, 1, :, :, 11], conditions[11], "G1")
    p22 = Eachdrug(tens[2, 1, :, :, 11], conditions[11], "S/G2")

    pl = plotGrid((6, 4), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, nothing, nothing];)
    return draw(SVG("figure1.svg", 18inch, 20inch), pl)
end

function figure2()
    setGadflyTheme()
    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)
    p1 = Eachdrug(tens[1, 1, :, :, 1] .+ tens[2, 1, :, :, 1], conditions[1], "total")
    p2 = Eachdrug(tens[1, 1, :, :, 2] .+ tens[2, 1, :, :, 2], conditions[2], "total")
    p3 = Eachdrug(tens[1, 1, :, :, 3] .+ tens[2, 1, :, :, 3], conditions[3], "total")
    p4 = Eachdrug(tens[1, 1, :, :, 4] .+ tens[2, 1, :, :, 4], conditions[4], "total")
    p5 = Eachdrug(tens[1, 1, :, :, 5] .+ tens[2, 1, :, :, 5], conditions[5], "total")
    p6 = Eachdrug(tens[1, 1, :, :, 6] .+ tens[2, 1, :, :, 6], conditions[6], "total")
    p7 = Eachdrug(tens[1, 1, :, :, 7] .+ tens[2, 1, :, :, 7], conditions[7], "total")
    p8 = Eachdrug(tens[1, 1, :, :, 8] .+ tens[2, 1, :, :, 8], conditions[8], "total")
    p9 = Eachdrug(tens[1, 1, :, :, 9] .+ tens[2, 1, :, :, 9], conditions[9], "total")
    p10 = Eachdrug(tens[1, 1, :, :, 10] .+ tens[2, 1, :, :, 10], conditions[10], "total")
    p11 = Eachdrug(tens[1, 1, :, :, 11] .+ tens[2, 1, :, :, 11], conditions[11], "total")

    pl = plotGrid((3, 4), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nothing];)
    return draw(SVG("figure2.svg", 18inch, 14inch), pl)
end
