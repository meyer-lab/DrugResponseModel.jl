""" figure5 plot single drug treatments from combination experiments. """

function plot_fig5(concs, g1, g1data, tite, G, subPlabel)
    time = LinRange(0.0, 95.0, 193)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = ["control" "$(concs[2]) nM" "$(concs[3]) nM" "$(concs[4]) nM" "$(concs[5]) nM"],
        fg_legend = :transparent,
        palette = :PuBu_5,
        title = tite,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        xlabel = "time [hr]",
        xticks = 0:24.0:96.0,
        ylabel = "$G cell number",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(time, g1data, lw = 2, linestyle = :dot, label = ["" "" "" "" ""])
    annotate!(-1.0, 2.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 2.5))
    p
end


function figure5()

    gt1, gt2 = DrugResponseModel.import_combination("AU01001");
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");

    GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
    GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

    meanGS1 = mean(GS1, dims = 4);
    meanGS2 = mean(GS2, dims = 4);
    meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims = 3), dims = 3)[:, :, 1]

    g1 = zeros(193, 5, 3)
    g2 = zeros(193, 5, 3)
    g1[:, 1, :] .= meanGS1[1, :, 1] # control
    g2[:, 1, :] .= meanGS1[2, :, 1]
    g1[:, 2:5, 1] .= meanGS1[1, :, 2:5]
    g2[:, 2:5, 1] .= meanGS1[2, :, 2:5]
    g1[:, 2:5, 2] .= meanGS1[1, :, 19:22]
    g2[:, 2:5, 2] .= meanGS1[2, :, 19:22]
    g1[:, 2:5, 3] .= meanGS1[1, :, 7:10]
    g2[:, 2:5, 3] .= meanGS1[2, :, 7:10]

    G1 = zeros(193, 5, 3)
    G2 = zeros(193, 5, 3)
    concs = zeros(5, 3)
    concs[:, [1, 3]] .= [0.0, 25.0, 50.0, 100.0, 250.0]
    concs[:, 2] .= [0.0, 5.0, 10.0, 17.0, 30.0]

    time = LinRange(0.0, 95.0, 193)
    p = parameters()
    ps = getODEparams(p, concs)
    for i=1:3
        for j=1:5
            G1[:, j, i], G2[:, j, i], _ = predict(ps[:, j, i], ps[:, 1, i], time)
        end
    end

    p1 = plot_fig5(concs[:, 1], G1[:, :, 1], g1[:, :, 1], "lapatinib", "G1", "a")
    p2 = plot_fig5(concs[:, 2], G1[:, :, 2], g1[:, :, 2], "gemcitabine", "G1", "b")
    p3 = plot_fig5(concs[:, 3], G1[:, :, 3], g1[:, :, 3], "palbociclib", "G1", "c")
    p4 = plot_fig5(concs[:, 1], G2[:, :, 1], g2[:, :, 1], "lapatinib", "G2", "d")
    p5 = plot_fig5(concs[:, 2], G2[:, :, 2], g2[:, :, 2], "gemcitabine", "G2", "e")
    p6 = plot_fig5(concs[:, 3], G2[:, :, 3], g2[:, :, 3], "palbociclib", "G2", "f")
    p = plot(p1, p2, p3, p4, p5, p6, size=(1000, 700), layout=(2, 3))
    savefig(p, "figure5.svg")
end
