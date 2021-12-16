""" Cluster drugs based on their inferred ODE effects. """
drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]
params = ["α_1", "α_2", "α_3", "α_4", "β_1", "β_2", "β_3", "β_4", "γ_11", "γ_12", "γ_13", "γ_14", "γ_21", "γ_22", "γ_23", "γ_24"]

function pca_()
    scores = []
    loadings = []
    for i=1:8
        ps = DrugResponseModel.out_ODEparams()[:, :, i]' # 16 x 11
        # mean-center each feature
        for i = 1:11
            ps[:, i] = ps[:, i] .- mean(ps[:, i])
            ps[:, i] = ps[:, i] ./ std(ps[:, i])
        end
        M = fit(PCA, ps; pratio=1, maxoutdim=4)
        push!(scores, MultivariateStats.transform(M, ps))
        push!(loadings, MultivariateStats.projection(M)')
    end
    return scores, loadings
end

function scores_or_loadings_plot(i, Cls)

    all_pt = pca_()[i]
    df = [DataFrame(x=ps_tr[1, :], y=ps_tr[2, :], class=Cls) for ps_tr in all_pt]
    plots = []
    for i=1:8
        push!(plots, Gadfly.plot(df[i], x="x", y="y", Geom.point, label="class", Geom.label, Guide.title("concentration $i"), Guide.xlabel("PC1"), Guide.ylabel("PC2")))
    end
    return plots
end

function figure3()
    setGadflyTheme()

    ppp = scores_or_loadings_plot(1, drugs)
    pl = plotGrid((2, 4), [ppp...];)
    return draw(SVG("figure3a.svg", 16inch, 8inch), pl)
end

function figure31()
    setGadflyTheme()

    ppp = scores_or_loadings_plot(2, params)
    pl = plotGrid((2, 4), [ppp...];)
    return draw(SVG("figure3b.svg", 12inch, 4inch), pl)
end