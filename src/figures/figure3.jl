""" Cluster drugs based on their inferred ODE effects. """
drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

function kmeans_()
    ps = DrugResponseModel.out_ODEparams()[:, :, 5]
    results = []
    for i = 2:6
        push!(results, kmeans(ps', i))
    end
    return results, ps
end

function clusterings(res)
    cls = []
    for i=1:11
        add_on = "c_"
        class = string(res.assignments[i])
        push!(cls, add_on * class)
    end
    df = DataFrame(x=drugs, y=res.assignments, class=cls)
    p = Gadfly.plot(df, x="x", y="y", Geom.point, color="class", Guide.xlabel("drugs"), Guide.ylabel("clusters"))
    return p
end

function pca_()
    ps2 = DrugResponseModel.out_ODEparams()[:, :, 5]
    M = fit(PCA, ps2'; pratio=1, maxoutdim=4)
    ps_transformed = MultivariateStats.transform(M, ps2')
    print(principalvars(M) ./ tvar(M) * 100)
return ps_transformed
end

function pca_plot()

    ps_tr = pca_()
    df = DataFrame(x=ps_tr[1, :], y=ps_tr[2, :], class=drugs)
    p = Gadfly.plot(df, x="x", y="y", Geom.point, label="class", Geom.label, Guide.xlabel("PC1"), Guide.ylabel("PC2"))
    return p
end

function figure3()
    setGadflyTheme()

    res, ps = kmeans_()
    plts = []
    for i=1:5
        push!(plts, clusterings(res[i]))
    end
    push!(plts, pca_plot())

    pl = plotGrid((2, 3), [plts...];)
    return draw(SVG("figure3.svg", 12inch, 8inch), pl)
end