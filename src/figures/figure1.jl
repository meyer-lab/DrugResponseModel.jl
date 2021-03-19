""" Figure 1 including model cartoon, time series simulations and parameters."""
# remember: to load the simple ODE params do: JLD.load("G1_simpleODE.jld")["data"]
function SSE(G1, G2, g1m, g2m, subPlabel)
    SSEs = zeros(5, 2)
    for i = 1:5
        SSEs[i, 2] = norm(G1[:, :, i] - g1m[:, 1:7, i]) + norm(G2[:, :, i] - g2m[:, 1:7, 1])
    end
    SSEs[:, 1] .= [379.145, 449.709, 483.859, 413.068, 402.336]
    ctg = repeat(["w/o LCT", "w/ LCT"], inner = 5)
    nam = repeat(["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"], outer = 2)

    groupedbar(
        nam,
        SSEs,
        group = ctg,
        xrotation = 30,
        xlabel = "Drugs",
        ylabel = "SSE",
        title = "Sum of Squared Errors",
        bar_width = 0.45,
        lw = 0,
        framestyle = :box,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.25cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.85cm,
        right_margin = 1.25cm,
    )
    annotate!(-1, 555.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
end

function plot_fig1(concs, g1, g1data, tite, G, subPlabel)
    time = LinRange(0.0, 95.0, 189)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = ["control" "$(concs[2]) nM" "$(concs[3]) nM" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM"],
        fg_legend = :transparent,
        palette = :PuBu_7,
        title = tite,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        xlabel = "time [hr]",
        ylabel = "$G cell number",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(time, g1data, lw = 2, linestyle = :dot, label = ["" "" "" "" "" "" ""])
    annotate!(-25.0, 57.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 50))
    p
end

function plot_pG1(efcs, ymax, drName, ylabel, subPlabel, plus)

    x = ["G11", "G12", "G21", "G22", "G23", "G24"]
    y1 = efcs[1:6, 1]
    y2 = efcs[1:6, 8]

    scatter(
        x,
        y1,
        color = "cyan4",
        xlabel = "sub-phase",
        label = "control",
        markerstrokewidth = 0,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$drName effects",
    )
    scatter!(
        x,
        y2,
        color = "cyan3",
        xlabel = "sub-phase",
        label = "max effect",
        markerstrokewidth = 0,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$drName effects",
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.05, ymax))
end

function plot_pG2(efcs, ymax, drName, ylabel, plus)
    x = ["G11", "G12", "", "G21", "G22", "G23", "G24"]
    y1 = [efcs[1:2, 1], nothing, efcs[3:6, 1]]
    y2 = [efcs[1:2, 8], nothing, efcs[3:6, 8]]

    scatter(
        x,
        y1,
        color = "cyan4",
        xlabel = "sub-phase",
        label = "control",
        markerstrokewidth = 0,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$drName effects",
    )
    scatter!(
        x,
        y2,
        color = "cyan3",
        xlabel = "sub-phase",
        label = "max effect",
        markerstrokewidth = 0,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$drName effects",
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.05, ymax))
end

function figure1()

    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2

    ps = [
        44.184,
        1.24076,
        0.0692788,
        0.0460918,
        0.3822,
        0.854034,
        0.605391,
        0.771326,
        0.0138293,
        0.00183699,
        0.000293753,
        0.0127534,
        0.00011816,
        0.0142541,
        60.6069,
        0.899573,
        1.99993,
        0.0748216,
        1.99326,
        0.468332,
        1.99864,
        1.22536,
        0.000141615,
        0.0318616,
        0.000216899,
        8.80158e-7,
        0.598489,
        0.00110572,
        6.68492,
        2.05974,
        1.99936,
        0.167588,
        0.507586,
        0.316074,
        0.248084,
        0.826596,
        1.6164e-5,
        3.10987e-6,
        3.55996e-5,
        7.73526e-6,
        0.0774056,
        8.26708e-5,
        3.34656,
        2.83739,
        0.0907361,
        0.108245,
        1.9758,
        1.96985,
        1.9993,
        0.210137,
        0.0690636,
        1.30442e-5,
        0.0767181,
        0.00991078,
        6.87891e-5,
        1.45086e-5,
        18.2253,
        1.1841,
        1.00505,
        0.0735852,
        1.97326,
        0.783828,
        0.45769,
        1.99355,
        0.0519941,
        0.000533671,
        0.00204743,
        9.52975e-5,
        5.23806e-5,
        0.0677505,
        0.339953,
        0.403341,
        0.802518,
        0.470576,
        1.298,
        0.423103,
    ]
    efcs = getODEparams(ps, concs)

    # ******* model simulations ********
    G1 = zeros(189, 7, 5)
    G2 = zeros(189, 7, 5)

    t = LinRange(0.0, 95.0, 189)
    for k = 1:5 # drug number
        for i = 1:7 # concentration number
            G1[:, i, k], G2[:, i, k], _ = predict(efcs[:, i, k], efcs[:, 1, k], t)
        end
    end

    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_fig1(concs[:, 1], G1[:, :, 1], g1m[:, 1:7, 1, 1], "Lapatinib", "G1", "b")
    p2 = plot_fig1(concs[:, 1], G2[:, :, 1], g2m[:, 1:7, 1, 1], "Lapatinib", "G2", "c")
    p3 = plot_fig1(concs[:, 2], G1[:, :, 2], g1m[:, 1:7, 2, 1], "Doxorubicin", "G1", "d")
    p4 = plot_fig1(concs[:, 2], G2[:, :, 2], g2m[:, 1:7, 2, 1], "Doxorubicin", "G2", "e")
    p5 = SSE(G1, G2, g1m, g2m, "f")
    p6 = plot_pG1(efcs[1:6, :, 1], 2.0, "Lapatinib", "progression rates", "g", 0.3)
    p7 = plot_pG1(efcs[7:12, :, 1], 0.1, "Lapatinib", "death rates", "h", 0.02)
    p8 = plot_pG1(efcs[1:6, :, 2], 2.0, "Doxorubicin", "progression rates", "i", 0.3)
    p9 = plot_pG1(efcs[7:12, :, 2], 0.1, "Doxorubicin", "death rates", "j", 0.02)

    figure1 = plot(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, size = (2000, 700), layout = (2, 5))
    #annotate!(-70, 30.0, text("a", :black, :left, Plots.font("Helvetica Bold", 15)))
    savefig(figure1, "figure1.svg")
end
