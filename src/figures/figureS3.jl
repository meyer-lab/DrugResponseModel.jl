"""To plot Supplementary Figure 3 A, the Erlang distributions/"""

function figure300()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2
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
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
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
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        bottom_margin = 1.5cm,
        alpha = 0.6,
        lw = 3,
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
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
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
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        bottom_margin = 1.5cm,
        lw = 3,
        alpha = 0.6,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )

    figs1 = plot(fig1, fig2, layout = (2, 1), size = (800, 400))
    savefig(figs1, "SupplementaryFigure2.svg")
end