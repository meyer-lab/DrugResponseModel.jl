""" Figure 1 including model cartoon, time series simulations and parameters."""
# remember: to load the simple ODE params do: JLD.load("G1_simpleODE.jld")["data"]
function SSE(G1, G2, g1m, g2m, subPlabel)
    SSEs = zeros(5, 2)
    G1ref = JLD.load("data/G1ref.jld")["G1ref"]
    G2ref = JLD.load("data/G2ref.jld")["G2ref"]
    for i = 1:5
        SSEs[i, 1] = norm(G1ref[:, :, i] - g1m[:, :, i]) + norm(G2ref[:, :, i] - g2m[:, :, 1])
        SSEs[i, 2] = norm(G1[:, :, i] - g1m[:, :, i]) + norm(G2[:, :, i] - g2m[:, :, 1])
    end
    ctg = repeat(["w/o LCT", "w/ LCT"], inner = 5)
    nam = repeat(["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"], outer = 2)

    groupedbar(
        nam,
        SSEs,
        group = ctg,
        xrotation = 30,
        xlabel = "drugs",
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

function plot_fig1(concs, g1, g1data, tite, G, subPlabel, palet)
    time = LinRange(0.0, 95.0, 189)

    p = plot(
        time,
        g1,
        lw = 2,
        legend = :topleft,
        label = ["control" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM" "$(concs[8]) nM"],
        fg_legend = :transparent,
        palette = palet,
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
    plot!(time, g1data, lw = 2, linestyle = :dot, label = ["" "" "" "" "" "" ""])
    annotate!(-1.0, 2.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 2.5))
    p
end

function plot_pG1_mean(y1, y2, ymax, Phasename, ylabel, subPlabel, plus)

    x = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]
    conts = scatter(
        x,
        y1',
        color = "red",
        xlabel = "drugs",
        xrotation = 30,
        label = "Control",
        markerstrokewidth = 0,
        markersize = 8,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$Phasename effects",
    )
    scatter!(
        x,
        y2',
        color = "cyan4",
        xlabel = "drugs",
        label = "Emax",
        markerstrokewidth = 0,
        markersize = 8,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
        title = "$Phasename effects",
    )
    annotate!(-0.5, (ymax + plus), text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.1, ymax))
end


function figure3()

    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2

    ps = parameters()
    efcs = getODEparams(ps, concs)

    # ******* model simulations ********
    G1 = zeros(189, 8, 5)
    G2 = zeros(189, 8, 5)

    t = LinRange(0.0, 95.0, 189)
    for k = 1:5 # drug number
        for i = 1:8 # concentration number
            G1[:, i, k], G2[:, i, k], _ = predict(efcs[:, i, k], efcs[:, 1, k], t)
        end
    end

    # deaths
    deathContG1 = zeros(2, 5)
    deathEmaxG1 = zeros(2, 5)
    deathContG1[1, :] .= (efcs[7, 1, :]) ./ (efcs[7, 1, :] .+ efcs[1, 1, :])
    deathContG1[2, :] .= (1 .- deathContG1[1, :]) .* (efcs[8, 1, :]) ./ (efcs[8, 1, :] .+ efcs[2, 1, :])
    deathEmaxG1[1, :] .= (efcs[7, 8, :]) ./ (efcs[7, 8, :] .+ efcs[1, 8, :])
    deathEmaxG1[2, :] .= (1 .- deathEmaxG1[1, :]) .* (efcs[8, 8, :]) ./ (efcs[8, 8, :] .+ efcs[2, 8, :])
    deathContG2 = zeros(4, 5)
    deathEmaxG2 = zeros(4, 5)
    deathContG2[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[3, 1, :])
    deathEmaxG2[1, :] .= (efcs[9, 8, :]) ./ (efcs[9, 8, :] .+ efcs[3, 8, :])
    for i = 10:12
        deathContG2[i - 8, :] = (1 .- deathContG2[i - 9, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 6, 1, :])
        deathEmaxG2[i - 9, :] = (1 .- deathEmaxG2[i - 9, :]) .* (efcs[i, 8, :]) ./ (efcs[i, 8, :] .+ efcs[i - 6, 8, :])
    end

    G1short = zeros(189, 6, 5)
    G2short = zeros(189, 6, 5)
    g1mshort = zeros(189, 6, 5)
    g2mshort = zeros(189, 6, 5)
    G1short[:, 1, :] .= G1[:, 1, :]
    G2short[:, 1, :] .= G2[:, 1, :]
    g1mshort[:, 1, :] .= g1m[:, 1, :]
    g2mshort[:, 1, :] .= g2m[:, 1, :]
    G1short[:, 2:6, :] .= G1[:, 4:8, :]
    G2short[:, 2:6, :] .= G2[:, 4:8, :]
    g1mshort[:, 2:6, :] .= g1m[:, 4:8, :]
    g2mshort[:, 2:6, :] .= g2m[:, 4:8, :]
    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_fig1(concs[:, 1], G1short[:, :, 1], g1mshort[:, :, 1, 1], "Dynamical Model Fits - Lapatinib", "G1", "B", :PuBu_6)
    p2 = plot_fig1(concs[:, 1], G2short[:, :, 1], g2mshort[:, :, 1, 1], "Dynamical Model Fits - Lapatinib", "S/G2", "C", :PuBu_6)
    p3 = plot_fig1(concs[:, 3], G1short[:, :, 3], g1mshort[:, :, 3, 1], "Dynamical Model Fits - Gemcitabine", "G1", "D", :PuBu_6)
    p4 = plot_fig1(concs[:, 3], G2short[:, :, 3], g2mshort[:, :, 3, 1], "Dynamical Model Fits - Gemcitabine", "S/G2", "E", :PuBu_6)
    p5 = plot_pG1_mean(mean(efcs[1:2, 1, :], dims = 1), mean(efcs[1:2, 8, :], dims = 1), 2.0, "G1 ", "progression rates [1/hr]", "F", 0.3)
    p6 = plot_pG1_mean(mean(efcs[3:6, 1, :], dims = 1), mean(efcs[3:6, 8, :], dims = 1), 2.0, "S/G2 ", "progression rates [1/hr]", "G", 0.3)
    p7 = plot_pG1_mean(sum(deathContG1, dims = 1), sum(deathEmaxG1, dims = 1), 1.0, "G1", "death probability", "H", 0.06)
    p8 = plot_pG1_mean(sum(deathContG2, dims = 1), sum(deathEmaxG2, dims = 1), 1.0, "S/G2", "death probability", "I", 0.06)
    p9 = SSE(G1[:, :, :], G2[:, :, :], g1m, g2m, "J")

    figure1 = plot(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, size = (2200, 700), layout = (2, 5))
    savefig(figure1, "figure3.svg")
end

# replicates' parameters:
# p1 = [59.132, 1.23216, 2.93162, 0.0224589, 0.275675, 2.99447, 0.493825, 0.582667, 0.0014566, 0.00738793, 0.0297349, 0.0763022, 0.00317598, 0.000262541, 172.905, 0.8691, 2.94725, 0.133166, 0.234383, 2.73361, 0.464481, 2.98659, 0.000172721, 0.000107943, 0.202668, 0.00391187, 0.000557085, 0.345508, 3.52234, 1.74014, 0.000323184, 0.255954, 0.100639, 0.124028, 0.730151, 0.901477, 0.0202751, 1.83514e-5, 1.5387e-5, 0.0315427, 0.000701986, 3.44425e-5, 2.37072, 2.10702, 0.542407, 0.0866461, 1.66253, 1.64308, 0.199229, 1.82007, 0.0594886, 0.0267128, 0.000201773, 0.000405221, 0.0151659, 0.00142696, 167.719, 0.684362, 2.18003, 0.0284909, 1.37066, 1.85981, 0.267127, 2.99785, 0.00214824, 3.36141e-5, 0.134529, 0.000727143, 6.67082e-5, 0.126732, 2.99503, 0.169147, 0.429907, 2.98509, 0.471941, 0.474202]
# p2 = [42.146, 1.25318, 0.00756712, 0.0697291, 0.085851, 0.00123844, 0.775785, 0.586675, 0.00613917, 0.000414688, 0.00326492, 7.08622e-5, 0.000116708, 0.0250105, 35.0601, 0.713745, 0.930779, 0.00027068, 1.7175, 1.9972, 0.223518, 1.99703, 0.00021207, 0.0645482, 0.0631564, 0.376948, 0.000373298, 0.000581707, 5.05173, 2.3473, 0.502171, 0.348143, 0.414099, 0.203361, 0.204054, 1.98124, 2.67144e-5, 1.9479e-5, 4.97218e-5, 0.0645767, 1.52479e-5, 7.28062e-5, 3.23095, 2.09279, 0.0897522, 0.0948307, 1.97429, 0.0507758, 0.921704, 0.0677673, 9.64556e-6, 0.188601, 0.00153864, 7.74546e-5, 3.10941e-6, 0.032831, 26.8248, 1.09054, 0.0957954, 0.792823, 0.787165, 1.27552, 0.356341, 1.9483, 1.98233e-5, 0.0238753, 0.000176911, 0.00200283, 0.021934, 0.00858807, 0.244592, 1.26283, 1.72032, 1.9972, 0.292368, 0.775993]
# p3 = [57.7701, 0.95457, 0.0325867, 0.0412113, 0.429698, 0.967497, 0.849089, 2.9987, 0.0149979, 5.57127e-6, 4.39279e-5, 2.01576e-5, 3.30843e-5, 0.0629629, 88.7115, 1.01634, 2.73977, 0.0700161, 2.48209, 0.564991, 2.99984, 2.99375, 3.49147e-5, 0.0289362, 8.14621e-5, 2.44595e-5, 0.70729, 0.000131428, 7.78495, 1.88809, 2.99927, 0.177006, 0.254106, 0.173401, 1.66882, 2.99713, 1.17434e-5, 7.03714e-6, 6.09986e-6, 0.0504083, 6.46098e-5, 9.38154e-5, 4.13404, 6.6356, 0.188563, 0.0646219, 2.99742, 0.482296, 0.286644, 0.0666906, 2.57614e-6, 2.99968e-6, 0.858529, 1.39496e-5, 4.66364e-6, 0.0409683, 26.1721, 1.12445, 1.88882, 0.0624009, 2.0323, 0.313346, 2.70712, 2.99926, 0.0586253, 0.000730662, 0.0784669, 4.3896e-5, 0.000300769, 0.06418, 0.319625, 0.431269, 0.386734, 0.618726, 0.630307, 1.01177]

function p_t()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2
    gt = g1m .+ g2m
    time = LinRange(0.0, 95.0, 189)

    p = plot(
        time,
        gt[:, end, 1],
        label = "lapatinib",
        xlabel = "time [hr]",
        ylabel = "cell numbers",
        fg_legend = :transparent,
        lw = 2,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    p = plot!(
        time,
        gt[:, end, 2],
        label = "doxorubicin",
        xlabel = "time [hr]",
        ylabel = "cell numbers",
        fg_legend = :transparent,
        lw = 2,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        time,
        gt[:, end, 3],
        label = "gemcitabine",
        xlabel = "time [hr]",
        ylabel = "cell numbers",
        fg_legend = :transparent,
        lw = 2,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        time,
        gt[:, end, 4],
        label = "paclitaxel",
        xlabel = "time [hr]",
        ylabel = "cell numbers",
        fg_legend = :transparent,
        lw = 2,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        xticks = 0:24.0:96.0,
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        time,
        gt[:, end, 5],
        label = "palbociclib",
        xlabel = "time [hr]",
        ylabel = "cell numbers",
        fg_legend = :transparent,
        lw = 2,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        xticks = 0:24.0:96.0,
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    ylims!((0.0, 40.0))
    savefig(p, "tot.svg")
end

function temp()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2
    gt = g1m .+ g2m
    time = LinRange(0.0, 95.0, 189)
    p1 = scatter(time, g1m[:, 6, 1], xlabel = "time [hr]", ylabel = "G1 cell counts")
    p2 = scatter(time, g2m[:, 6, 1], xlabel = "time [hr]", ylabel = "S/G2 cell counts")
    pp = plot(p1, p2, size = (700, 300), layout = (1, 2))
    savefig(pp, "temp.svg")
end
