""" To plot the Erlang distributions. """

function plotErlang()
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    path_g = string("/", basePath, "/CCDurations.xlsx")
    columns, labels = XLSX.readtable(path_g, "CONTROL")
    columns = hcat(columns)

    # G1 fit
    G1 = columns[:, 1][1]
    G1 = convert(Array{Float64, 1}, G1)
    # Gamma{Float64}(α=6.6035911662388775, θ=2.7246247139426316)
    d1 = Gamma(6.6, 2.72)
    lo, hi = quantile.(d1, [0.0001, 0.9999])
    x1 = range(lo, hi; length = length(G1))
    pdfG1 = pdf.(d1, x1)
    coeffG1 = maximum(G1) / maximum(pdfG1)

    # G2 fit
    Distributions.fit_mle(Gamma, G1)
    G2 = columns[:, 1][2][1:514]
    G2 = convert(Array{Float64, 1}, G2)
    Distributions.fit_mle(Gamma, G2)
    # Gamma{Float64}(α=22.59370017485455, θ=0.9964125282786628)
    d2 = Gamma(22.6, 1.0)
    lo, hi = quantile.(d2, [0.0001, 0.9999])
    x2 = range(lo, hi; length = length(G2))
    pdfG2 = pdf.(d2, x2)
    coeffG2 = maximum(G2) / maximum(pdfG2)
    # plots
    fig1 = histogram(
        G1,
        normalize = false,
        xlabel = "time[hr]",
        ylabel = "abundance",
        label = "Control data",
        title = "G1 phase lengths",
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        alpha = 0.6,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        x1,
        coeffG1 .* pdfG1,
        label = "Gamma dist.",
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        alpha = 0.6,
        lw = 2,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    fig2 = histogram(
        G2,
        normalize = false,
        label = "Control data",
        xlabel = "time[hr]",
        ylabel = "abundance",
        title = "S/G2 phase lengths",
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        alpha = 0.6,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        x2,
        coeffG2 .* pdfG2,
        label = "Gamma dist.",
        density = false,
        titlefont = Plots.font("Helvetica", 12),
        legendfont = Plots.font("Helvetica", 9),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        lw = 2,
        alpha = 0.6,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    figs2 = plot(fig1, fig2, size = (900, 350))
    savefig(figs2, "figureS2.svg")
end

function output_durations()
    concs, _, _, _ = load(189, 1)
    ps = [
        51.0122,
        1.19478,
        0.0123853,
        0.197453,
        0.783039,
        6.53136e-5,
        1.35692e-6,
        0.284673,
        0.00521293,
        3.69958e-7,
        0.00913979,
        0.0258875,
        3.04229e-6,
        0.00527735,
        18.4107,
        1.38004,
        0.288625,
        9.6902e-9,
        0.787761,
        1.02151,
        1.99999,
        0.106618,
        4.35605e-9,
        0.0478454,
        1.22383e-7,
        1.04499e-7,
        0.381662,
        2.39835e-9,
        4.75582,
        1.78552,
        0.481014,
        0.404215,
        0.471125,
        0.187735,
        1.99999,
        0.255864,
        1.35294e-9,
        7.07919e-9,
        1.74332e-9,
        0.0672485,
        4.87662e-8,
        4.45473e-9,
        7.0734,
        2.47932,
        0.066145,
        5.62597e-8,
        1.94036,
        2.0,
        2.0,
        0.00866935,
        1.22435e-9,
        9.23547e-7,
        2.0,
        2.14921e-7,
        1.23361e-7,
        0.0174862,
        36.8515,
        1.11516,
        0.0806277,
        0.726529,
        1.92473,
        1.99999,
        1.97768,
        0.319934,
        2.65382e-9,
        6.12668e-9,
        0.0197645,
        1.06389e-6,
        5.28303e-8,
        0.0308013,
        0.196915,
        2.0,
        1.92313,
        2.0,
        1.99921,
        0.199044,
    ]
    efcs = getODEparams(ps, concs)

    gi = zeros(2, 8, 5)
    gi[1, :, :] .= (4 ./ efcs[1, :, :] .+ 4 ./ efcs[2, :, :])
    gi[2, :, :] .= (5 ./ efcs[3, :, :] .+ 5 ./ efcs[4, :, :] .+ 5 ./ efcs[5, :, :] .+ 5 ./ efcs[6, :, :])

    df1 = DataFrames.DataFrame(lap = gi[1, :, 1], dox = gi[1, :, 2], gem = gi[1, :, 3], pax = gi[1, :, 4], palbo = gi[1, :, 5])
    XLSX.writetable("durationsG1.xlsx", df1)
    df2 = DataFrames.DataFrame(lap = gi[2, :, 1], dox = gi[2, :, 2], gem = gi[2, :, 3], pax = gi[2, :, 4], palbo = gi[2, :, 5])
    XLSX.writetable("durationsG2.xlsx", df2)
end
