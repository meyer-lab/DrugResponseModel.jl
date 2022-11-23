""" Supplementary figures for paper, new cell lines. """

function plot_fig1_data(concs, g1data, tite, G, subPlabel, palet, time)
    p = Plots.plot(
        time,
        g1data,
        lw = 5,
        legend = :topleft,
        label = ["control" "$(concs[2]) nM" "$(concs[3]) nM" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM" "$(concs[8]) nM"],
        fg_legend = :transparent,
        palette = palet,
        title = tite,
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        xlabel = "Time [hr]",
        xticks = 0:24.0:96.0,
        ylabel = "Normalized Cell Counts",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    annotate!(-0.5, 1.5, Plots.text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 4))
    p
end

function plot_fig1_percG1(concs, g1data, tite, G, subPlabel, palet, time)
    p = Plots.plot(
        time,
        g1data,
        lw = 5,
        legend = :topleft,
        label = ["control" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM" "$(concs[8]) nM"],
        fg_legend = :transparent,
        palette = palet,
        title = tite,
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        xlabel = "Time [hr]",
        xticks = 0:24.0:96.0,
        ylabel = "G1 Percentage",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    annotate!(-0.5, 1.5, Plots.text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 1.0))
    p
end

""" To plot the fits and accumulated cell death for each cell line, we do the following:
1. tensor, names, concs, conds = DrugResponseModel.__cellLineName__()
where __cellLineName__ could be one of [hcc_all, mt1_all, mda_all]
2. imporing the estimated parameters according to the cell line, one of [ps_hcc, ps_mt1, ps_mda] above.
3. DrugResponseModel.figure70(tensor, names, concs, conds, ps)"""

function figureS1()
    ENV["GKSwstype"]="nul"
    tensor, names, concs, conds = DrugResponseModel.mda_all()
    t = LinRange(0.0, 96, size(tensor)[2])

    gmshort = zeros(size(tensor)[2], 6, 6, 2) # datapoints x concs x drugs x g1/g2
    for i=1:2
        gmshort[:, 1, :, i] .= tensor[i, :, 1, :]
        gmshort[:, 2:, :, i] .= tensor[i, 4:, :, :]
    end


    p1 = plot_fig1_data(concs[1], gmshort[:, :, 1, 1] .+ gmshort[:, :, 1, 2], string("MDA-MB-157 treated with ", names[1]), "", "A", :PuBu_8, t)
    p2 = plot_fig1_data(concs[2], gmshort[:, :, 2, 1] .+ gmshort[:, :, 2, 2], string("MDA-MB-157 treated with ", names[2]), "", "", :PuBu_8, t)
    p3 = plot_fig1_data(concs[3], gmshort[:, :, 3, 1] .+ gmshort[:, :, 3, 2], string("MDA-MB-157 treated with ", names[3]), "", "", :PuBu_8, t)
    p4 = plot_fig1_data(concs[4], gmshort[:, :, 4, 1] .+ gmshort[:, :, 4, 2], string("MDA-MB-157 treated with ", names[4]), "", "", :PuBu_8, t)
    p5 = plot_fig1_data(concs[5], gmshort[:, :, 5, 1] .+ gmshort[:, :, 5, 2], string("MDA-MB-157 treated with ", names[5]), "", "", :PuBu_8, t)
    p6 = plot_fig1_data(concs[6], gmshort[:, :, 6, 1] .+ gmshort[:, :, 6, 2], string("MDA-MB-157 treated with ", names[6]), "", "", :PuBu_8, t)

    p7 = plot_fig1_percG1(concs[1], gmshort[:, :, 1, 1] ./ (gmshort[:, :, 1, 1] .+ gmshort[:, :, 1, 2]), string("MDA-MB-157 treated with ", names[1]), "", "B", :PuBu_8, t)
    p8 = plot_fig1_percG1(concs[2], gmshort[:, :, 2, 1] ./ (gmshort[:, :, 2, 1] .+ gmshort[:, :, 2, 2]), string("MDA-MB-157 treated with ", names[2]), "", "", :PuBu_8, t)
    p9 = plot_fig1_percG1(concs[3], gmshort[:, :, 3, 1] ./ (gmshort[:, :, 3, 1] .+ gmshort[:, :, 3, 2]), string("MDA-MB-157 treated with ", names[3]), "", "", :PuBu_8, t)
    p10 = plot_fig1_percG1(concs[4], gmshort[:, :, 4, 1] ./ (gmshort[:, :, 4, 1] .+ gmshort[:, :, 4, 2]), string("MDA-MB-157 treated with ", names[4]), "", "", :PuBu_8, t)
    p11 = plot_fig1_percG1(concs[5], gmshort[:, :, 5, 1] ./ (gmshort[:, :, 5, 1] .+ gmshort[:, :, 5, 2]), string("MDA-MB-157 treated with ", names[5]), "", "", :PuBu_8, t)
    p12 = plot_fig1_percG1(concs[6], gmshort[:, :, 5, 1] ./ (gmshort[:, :, 6, 1] .+ gmshort[:, :, 6, 2]), string("MDA-MB-157 treated with ", names[6]), "", "", :PuBu_8, t)

    figure1 = Plots.plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (2600, 800), layout = (2, 6))
    Plots.savefig(figure1, "figureS1_MDAMB157.svg")
end
