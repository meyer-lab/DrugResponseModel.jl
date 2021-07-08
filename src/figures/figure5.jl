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
    concs = zeros(5, 3)
    concs[:, [1, 3]] .= [0.0, 25.0, 50.0, 100.0, 250.0]
    concs[:, 2] .= [0.0, 5.0, 10.0, 17.0, 30.0]

    G1 = zeros(193, 5, 3)
    G2 = zeros(193, 5, 3)

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

function plot_fig5s(concs, g1, g1data, tite, G, labels, subPlabel)
    time = LinRange(0.0, 95.0, 193)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = labels,
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

function combin()

    concs = zeros(5, 3)
    concs[:, [1, 3]] .= [0.0, 25.0, 50.0, 100.0, 250.0]
    concs[:, 2] .= [0.0, 5.0, 10.0, 17.0, 30.0]
    p = parameters()
    ps = getODEparams(p, concs)
    time = LinRange(0.0, 95.0, 193)

    lap_and_gem = DrugResponseModel.AllBliss_params(ps[:, :, 1], ps[:, :, 2]; n = 5)
    lap_and_palb = DrugResponseModel.AllBliss_params(ps[:, :, 1], ps[:, :, 3]; n = 5)
    gem_and_palb = DrugResponseModel.AllBliss_params(ps[:, :, 2], ps[:, :, 3]; n = 5)

    lap100_gems = zeros(3, 193, 5)
    lap100_plbs = zeros(3, 193, 5)
    palb50_laps = zeros(3, 193, 5)
    palb50_gems = zeros(3, 193, 5)
    gem10_laps = zeros(3, 193, 5)
    gem10_plbs = zeros(3, 193, 5)
    for i=1:5
        lap100_gems[1, :, i], lap100_gems[2, :, i], _ = predict(lap_and_gem[:, 3, i], lap_and_gem[:, 1, 1], time)
        lap100_plbs[1, :, i], lap100_plbs[2, :, i], _ = predict(lap_and_palb[:, 3, i], lap_and_palb[:, 1, 1], time)
        palb50_laps[1, :, i], palb50_laps[2, :, i], _ = predict(lap_and_palb[:, i, 2], lap_and_palb[:, 1, 1], time)
        palb50_gems[1, :, i], palb50_gems[2, :, i], _ = predict(gem_and_palb[:, i, 2], gem_and_palb[:, 1, 1], time)
        gem10_laps[1, :, i], gem10_laps[2, :, i], _ = predict(lap_and_gem[:, i, 2], lap_and_gem[:, 1, 1], time)
        gem10_plbs[1, :, i], gem10_plbs[2, :, i], _ = predict(gem_and_palb[:, 2, i], gem_and_palb[:, 1, 1], time)
    end

    lap100_gems[3, :, :] .= lap100_gems[1, :, :] .+ lap100_gems[2, :, :]
    lap100_plbs[3, :, :] .= lap100_plbs[1, :, :] .+ lap100_plbs[2, :, :]
    palb50_laps[3, :, :] .= palb50_laps[1, :, :] .+ palb50_laps[2, :, :]
    palb50_gems[3, :, :] .= palb50_gems[1, :, :] .+ palb50_gems[2, :, :]
    gem10_laps[3, :, :] .= gem10_laps[1, :, :] .+ gem10_laps[2, :, :]
    gem10_plbs[3, :, :] .= gem10_plbs[1, :, :] .+ gem10_plbs[2, :, :]

    gt1, gt2 = DrugResponseModel.import_combination("AU01001");
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");

    GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
    GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

    meanGS1 = mean(GS1, dims = 4);
    meanGS2 = mean(GS2, dims = 4);
    meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims = 3), dims = 3)[:, :, 1]

    g = zeros(3, 193, 6, 6)
    g[:, :, 1, :] .= meanGS2[:, :, 1, 1]
    ########### 1. Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
    g[:, :, 2, 1] .= meanGS1[:, :, 8, 1]
    g[:, :, 3:6, 1] .= meanGS2[:, :, 3:6, 1]
    # well 2: 3,4,5,6
    # palbo50 alone: well 1: 8
    
    ########### Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
    g[:, :, 2, 2] .= meanGS1[:, :, 8, 1]
    g[:, :, 3:6, 2] .= meanGS2[:, :, 21:24, 1]
    # well 2: 21, 22, 23, 24
    
    ########### Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM]
    g[:, :, 2, 3] .= meanGS1[:, :, 20, 1]
    g[:, :, 3, 3] .= meanGS1[:, :, 18, 1]
    g[:, :, 4, 3] .= meanGS1[:, :, 22, 1]
    g[:, :, 5:6, 3] .= meanGS2[:, :, 23:24, 1]
    # well 1: 18, 23, 24
    # gem 10 alone: well 1: 20
    
    ########### Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
    g[:, :, 2, 4] .= meanGS1[:, :, 20, 1]
    g[:, :, 3:6, 4] .= meanGS2[:, :, 9:12, 1]
    # well 2: 9, 10, 11, 12
    # gem 10 alone: well 1: 20
    
    
    ########### Lap 100 nM + gemcitabines [5, 10, 17 nM, 30 nM]
    g[:, :, 2, 5] .= meanGS1[:, :, 4, 1]
    g[:, :, 3, 5] .= meanGS2[:, :, 7, 1]
    g[:, :, 4, 5] .= meanGS2[:, :, 11, 1]
    g[:, :, 5, 5] .= meanGS2[:, :, 13, 1]
    g[:, :, 6, 5] .= meanGS2[:, :, 19, 1]
    # well 2: 7, 11, 13, 19
    # lap 100 alone well 1 : 4
    
    ########### Lap 100 nM + palbociclibs [25 nM, 50, 100 nM, 250 nM]
    g[:, :, 2, 6] .= meanGS1[:, :, 4, 1]
    g[:, :, 3, 6] .= meanGS2[:, :, 8, 1]
    g[:, :, 4, 6] .= meanGS2[:, :, 5, 1]
    g[:, :, 5, 6] .= meanGS2[:, :, 14, 1]
    g[:, :, 6, 6] .= meanGS2[:, :, 20, 1]

    p1 = plot_fig5s(concs[:, 2], lap100_gems[3, :, :], g[3, :, 2:6, 5], "lpt 100 + gemc", "total", ["lpt100" "lpt100+gem5" "lpt100+gem10" "lpt100+gem17" "lpt100+gem30"], "A")
    p2 = plot_fig5s(concs[:, 3], lap100_plbs[3, :, :], g[3, :, 2:6, 6], "lpt 100 + palbos", "total", ["lpt100" "lpt100+plb25" "lpt100+plb50" "lpt100+plb100" "lpt100+plb250"], "B")
    p3 = plot_fig5s(concs[:, 3], palb50_laps[3, :, :], g[3, :, 2:6, 1], "palbo 50 + lapts", "total", ["plb50" "plb50+lpt25" "plb50+lpt50" "plb50+lpt100" "plb50+lpt250"], "C")
    p4 = plot_fig5s(concs[:, 3], palb50_gems[3, :, :], g[3, :, 2:6, 2], "palbo 50 + gems", "total", ["plb50" "plb50+gem5" "plb50+gem10" "plb50+gem17" "plb50+gem30"], "D")
    p5 = plot_fig5s(concs[:, 3], gem10_laps[3, :, :], g[3, :, 2:6, 4], "gem10 + laps", "total", ["gem10" "gem10+lpt25" "gem10+lpt50" "gem10+lpt100" "gem10+lpt250"], "E")
    p6 = plot_fig5s(concs[:, 3], gem10_plbs[3, :, :], g[3, :, 2:6, 3], "gem10 + plbs", "total", ["gem10" "gem10+pl25" "gem10+plb50" "gem10+plb100" "gem10+plb250"], "F")

    pp = plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1000, 700))
    savefig(pp, "combin.svg")
end