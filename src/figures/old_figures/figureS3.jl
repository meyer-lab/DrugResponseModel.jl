""" To plot the Erlang distributions. """

function plotErlang()
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
    ps = parameters()
    efcs = getODEparams(ps, concs)

    # replace the second dimension of efcs with the ec50 effects
    for i = 1:5
        efcs[:, 2, i] = DrugResponseModel.EC50_params(ps, i)
    end

    G1 = zeros(189, 8, 5)
    G2 = zeros(189, 8, 5)

    t = LinRange(0.0, 95.0, 189)
    for k = 1:5 # drug number
        for i = 1:8 # concentration number
            G1[:, i, k], G2[:, i, k], _ = predict(efcs[:, i, k], efcs[:, 1, k], t)
        end
    end

    G1ref = JLD.load("data/G1ref.jld")["G1ref"]
    G2ref = JLD.load("data/G2ref.jld")["G2ref"]

    G1short = zeros(189, 6, 5)
    G2short = zeros(189, 6, 5)
    G1refshort = zeros(189, 6, 5)
    G2refshort = zeros(189, 6, 5)
    g1mshort = zeros(189, 6, 5)
    g2mshort = zeros(189, 6, 5)
    G1short[:, 1, :] .= G1[:, 1, :]
    G2short[:, 1, :] .= G2[:, 1, :]
    G1refshort[:, 1, :] .= G1ref[:, 1, :]
    G2refshort[:, 1, :] .= G2ref[:, 1, :]
    g1mshort[:, 1, :] .= g1m[:, 1, :]
    g2mshort[:, 1, :] .= g2m[:, 1, :]
    G1short[:, 2:6, :] .= G1[:, 4:8, :]
    G2short[:, 2:6, :] .= G2[:, 4:8, :]
    G1refshort[:, 2:6, :] .= G1ref[:, 3:7, :]
    G2refshort[:, 2:6, :] .= G2ref[:, 3:7, :]
    g1mshort[:, 2:6, :] .= g1m[:, 3:7, :]
    g2mshort[:, 2:6, :] .= g2m[:, 3:7, :]

    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p3 = DrugResponseModel.plot_fig1(concs[:, 1], G1refshort[:, :, 1], g1mshort[:, :, 1, 1], "Expon Model Fits - Lapatinib", "G1", "B", :YlOrBr_6)
    p4 = DrugResponseModel.plot_fig1(concs[:, 1], G2refshort[:, :, 1], g2mshort[:, :, 1, 1], "Expon Model Fits - Lapatinib", "S/G2", "C", :YlOrBr_6)
    p5 = DrugResponseModel.plot_fig1(concs[:, 3], G1refshort[:, :, 3], g1mshort[:, :, 3, 1], "Expon Model Fits - Gemcitabine", "G1", "D", :YlOrBr_6)
    p6 = DrugResponseModel.plot_fig1(concs[:, 3], G2refshort[:, :, 3], g2mshort[:, :, 3, 1], "Expon Model Fits - Gemcitabine", "S/G2", "E", :YlOrBr_6)
    figs2 = plot(p0, p0, p3, p4, p5, p6, fig1, fig2, layout = (4, 2), size = (800, 1700))
    savefig(figs2, "figureS3.svg")
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
