""" To plot the Erlang distributions. """

function plotErlang()
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    path_g = string("/", basePath, "/CCDurations.xlsx")
    df = DataFrame(XLSX.readtable(path_g, "CONTROL"))

    # G1 fit
    G1 = df[:, 1][1]
    G1 = convert(Array{Float64, 1}, G1)
    # Gamma{Float64}(α=6.6035911662388775, θ=2.7246247139426316)
    d1 = Gamma(6.6, 2.72)
    lo, hi = quantile.(d1, [0.0001, 0.9999])
    x1 = range(lo, hi; length = length(G1))
    pdfG1 = pdf.(d1, x1)
    coeffG1 = maximum(G1) / maximum(pdfG1)

    # G2 fit
    Distributions.fit_mle(Gamma, G1)
    G2 = df[:, 1][2][1:514]
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
