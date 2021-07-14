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
    ps = parameters()
    efcs = getODEparams(ps, concs)

    gi = zeros(2, 8, 5)
    gi[1, :, :] .= (2 ./ efcs[1, :, :] .+ 2 ./ efcs[2, :, :] .+ 2 ./ efcs[3, :, :] .+ 2 ./ efcs[4, :, :])
    gi[2, :, :] .= (5 ./ efcs[5, :, :] .+ 5 ./ efcs[6, :, :] .+ 5 ./ efcs[7, :, :] .+ 5 ./ efcs[8, :, :])

    df1 = DataFrames.DataFrame(
        lapG1 = gi[1, :, 1],
        doxG1 = gi[1, :, 2],
        gemG1 = gi[1, :, 3],
        paxG1 = gi[1, :, 4],
        palboG1 = gi[1, :, 5],
        lapG2 = gi[2, :, 1],
        doxG2 = gi[2, :, 2],
        gemG2 = gi[2, :, 3],
        paxG2 = gi[2, :, 4],
        palboG2 = gi[2, :, 5],
    )
    XLSX.writetable("durations.xlsx", df1)
end

# d1 = DataFrames.DataFrame(controlG1=G1[:, 1, 1], lpt5_G1=G1[:, 2, 1], lpt10_G1=G1[:, 3, 1], lpt25_G1=G1[:, 4, 1], lpt50_G1=G1[:, 5, 1], lpt100_G1=G1[:, 6, 1], lpt250_G1=G1[:, 7, 1], lpt500_G1=G1[:, 8, 1], 
#                           controlG2=G2[:, 1, 1], lpt5_G2=G2[:, 2, 1], lpt10_G2=G2[:, 3, 1], lpt25_G2=G2[:, 4, 1], lpt50_G2=G2[:, 5, 1], lpt100_G2=G2[:, 6, 1], lpt250_G2=G2[:, 7, 1], lpt500_G2=G2[:, 8, 1], 
#                           control=G1[:, 1, 1].+G2[:, 1, 1], lpt5=G1[:, 2, 1].+G2[:, 2, 1], lpt10=G1[:, 3, 1].+G2[:, 3, 1], lpt25=G1[:, 4, 1].+G2[:, 4, 1], lpt50=G1[:, 5, 1].+G2[:, 5, 1], lpt100=G1[:, 6, 1].+G2[:, 6, 1], lpt250=G1[:, 7, 1].+G2[:, 7, 1], lpt500=G1[:, 8, 1].+G2[:, 8, 1])