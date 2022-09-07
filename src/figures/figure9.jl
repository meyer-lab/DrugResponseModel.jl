""" Only plot the controls from the HCC1143 cells. """

function figure9()
    _, _, _, _, cl1, cl2 = DrugResponseModel.hcc_all()

    t = LinRange(0.0, 95, 182)
    dfs = [DataFrame(x=t, y1=cl1[1, :, i], y2=cl2[1, :, i], y12=cl1[2, :, i], y22=cl2[2, :, i]) for i = 1:size(cl2)[3]]
    DF = vcat(dfs...)
    a = [Gadfly.plot(DF, x="x", y="y1", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("G1 exp1"), Coord.Cartesian(ymin=-0.05,ymax=2.0))]
    push!(a, Gadfly.plot(DF, x="x", y="y2", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("G1 exp2"), Coord.Cartesian(ymin=-0.05,ymax=2.0)))
    push!(a, Gadfly.plot(DF, x="x", y="y12", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("SG2 exp1"), Coord.Cartesian(ymin=-0.05,ymax=2.0)))
    push!(a, Gadfly.plot(DF, x="x", y="y22", Geom.line, Guide.xlabel("time [hr]"), Guide.ylabel("SG2 exp2"), Coord.Cartesian(ymin=-0.05,ymax=2.0)))
    pl = plotGrid((2, 2), [a...];)
    return draw(SVG("figure9.svg", 8inch, 8inch), pl)

end