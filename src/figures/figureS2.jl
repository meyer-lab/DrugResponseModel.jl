""" To plot the Erlang distributions. """

function plotErlang()
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    path_g = string("/", basePath, "/CCDurations.xlsx")
    df = DataFrame(XLSX.readtable(path_g, "CONTROL"))

    G1 = df[:, 1][1]
    G1 = convert(Array{Float64,1}, G1)
    # Gamma{Float64}(α=6.6035911662388775, θ=2.7246247139426316)
    g1_rand = rand(Gamma(7.0, 2.72), length(G1))

    Distributions.fit_mle(Gamma, G1)
    G2 = df[:, 1][2][1:514]
    G2 = convert(Array{Float64,1}, G2)
    Distributions.fit_mle(Gamma, G2)
    # Gamma{Float64}(α=22.59370017485455, θ=0.9964125282786628)
    g2_rand = rand(Gamma(20.0, 1.0), length(G2))

    fig1 = histogram(G1, normalize=false, xlabel="time[hr]", ylabel="abundance", label="Control data", title="G1 phase lengths", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.5cm,
    alpha=0.6,
    fg_legend = :transparent,
    top_margin = 1.5cm,
    left_margin = 1.25cm,
    right_margin = 1.25cm)
    histogram!(g1_rand, normalize=false, label="Gamma dist.", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.5cm,
    alpha=0.6,
    fg_legend = :transparent,
    top_margin = 1.5cm,
    left_margin = 1.25cm,
    right_margin = 1.25cm)
    fig2 = histogram(G2, normalize=false, label="Control data", xlabel="time[hr]", ylabel="abundance", title="S/G2 phase lengths", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.5cm,
    alpha=0.6,
    fg_legend = :transparent,
    top_margin = 1.5cm,
    left_margin = 1.25cm,
    right_margin = 1.25cm)
    histogram!(g2_rand, label="Gamma dist.", normalize=false, titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.5cm,
    alpha=0.6,
    fg_legend = :transparent,
    top_margin = 1.5cm,
    left_margin = 1.25cm,
    right_margin = 1.25cm)
    figs2 = plot(fig1, fig2, size=(900, 350))
    savefig(figs2, "figureS2.svg")
end
