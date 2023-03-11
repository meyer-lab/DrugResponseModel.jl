""" Supplementary figure 1 for paper. """

"""This is to plot the replicates of AU565 cells treated with first round of drugs, including Lapatinib, Doxorubicin, Gemcitabine, Paclitaxel, and Palbociclib."""
function figure100()
    time = LinRange(0.0, 95.0, 189)

    concs, _, g1s1, g2s1 = load(189, 1);
    _, _, g1s2, g2s2 = load(189, 2);
    _, _, g1s3, g2s3 = load(189, 3);

    drs = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]
    css = [["control" "25nM" "50nM" "100nM" "250nM" "500nM"],
    ["control" "20nM" "50nM" "125nM" "250nM" "500nM"],
    ["control" "2.5nM" "5nM" "10nM" "30nM" "100nM"],
    ["control" "2nM" "3nM" "5nM" "7.5nM" "15nM"],
    ["control" "25nM" "50nM" "100nM" "250nM" "500nM"]]

    total_ps = []
    for i = 1:5
    p1 = Plots.plot(
        time,
        g1s1[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep1 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "G1 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuRd_6,
    )
    ylims!((0, 2))

    p2 = Plots.plot(
        time,
        g2s1[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep1 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "SG2 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuBu_6,
    )
    ylims!((0, 2))

    p3 = Plots.plot(
        time,
        g1s2[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep2 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "G1 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuRd_6,
    )
    ylims!((0, 2))

    p4 = Plots.plot(
        time,
        g2s2[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep2 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "SG2 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuBu_6,
    )
    ylims!((0, 2))

    p5 = Plots.plot(
        time,
        g1s3[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep3 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "G1 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuRd_6,
    )
    ylims!((0, 2))

    p6 = Plots.plot(
        time,
        g2s3[:, [1, 4, 5, 6, 7, 8], i],
        title = string("rep3 ", drs[i]),
        titlefontsize = 8,
        markersize = 1.0,
        markerstrokewidth = 0.0,
        legend = :topleft,
        label = css[i],
        xlabel = "time [hr]",
        ylabel = "SG2 cell proportions",
        xguidefontsize = 8,
        yguidefontsize = 8,
        legendfontsize = 8,
        lw = 3, 
        palette = :PuBu_6,
    )
    ylims!((0, 2))
    pp1 = [p1, p2, p3, p4, p5, p6]
    push!(total_ps, pp1)
    end
    s = vcat(total_ps...)
    pp = Plots.plot(s..., size = (2400, 1800), layout = (5, 6))
    savefig(pp, "SupplementaryFigure1.svg")
end
