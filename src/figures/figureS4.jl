"""Functions to plot the Supplementary Figure 4."""


function plot_pG1(efcs, ymax, Phasename, ylabel, subPlabel, plus)

    x = ["G11", "G12", "G13", "G14", "G21", "G22", "G23", "G24"]
    y1 = efcs[1:8, 1]
    y2 = efcs[1:8, 2]
    scatter(
        x,
        y1,
        color = "red",
        xlabel = "sub-phase",
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
        y2,
        color = "cyan4",
        xlabel = "sub-phase",
        label = "E_EC50",
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

function figureS2()

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

    ec50 = zeros(16, 5)
    conc_ec50 = zeros((1, 5))
    conc_ec50[1, 1:5] = [57.3 11.2 8.95 2.4 30.1]
    ec50 = getODEparams(ps, conc_ec50)[:, 1, :]

    # replace the second dimension of efcs with the ec50 effects

    efcs[:, 2, :] = ec50

    duration_efcs = zeros(8, 2, 5)
    death_efcs = zeros(8, 2, 5)

    duration_efcs[1:4, :, :] .= 2 ./ efcs[1:4, 1:2, :]
    duration_efcs[5:8, :, :] .= 5 ./ efcs[5:8, 1:2, :]

    deathContG1 = zeros(4, 5)
    deathEC50G1 = zeros(4, 5)
    deathContG1[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[1, 1, :])
    deathEC50G1[1, :] .= (ec50[9, :]) ./ (ec50[9, :] .+ ec50[1, :])
    for i = 2:4
        deathContG1[i, :] .= (1 .- deathContG1[i - 1, :]) .* (efcs[i + 8, 1, :]) ./ (efcs[i + 8, 1, :] .+ efcs[i, 1, :])
        deathEC50G1[i, :] .= (1 .- deathEC50G1[i - 1, :]) .* (ec50[i + 8, :]) ./ (ec50[i + 8, :] .+ ec50[i, :])
    end
    deathContG2 = zeros(4, 5)
    deathEC50G2 = zeros(4, 5)
    deathContG2[1, :] .= (efcs[13, 1, :]) ./ (efcs[13, 1, :] .+ efcs[5, 1, :])
    deathEC50G2[1, :] .= (ec50[13, :]) ./ (ec50[13, :] .+ ec50[5, :])
    for i = 14:16
        deathContG2[i - 12, :] = (1 .- deathContG2[i - 13, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 8, 1, :])
        deathEC50G2[i - 12, :] = (1 .- deathEC50G2[i - 13, :]) .* (ec50[i, :]) ./ (ec50[i, :] .+ ec50[i - 8, :])
    end

    death_efcs[1:4, 1, :] .= deathContG1
    death_efcs[5:8, 1, :] .= deathContG2
    death_efcs[1:4, 2, :] .= deathEC50G1
    death_efcs[5:8, 2, :] .= deathEC50G2
    # ******* model simulations ********
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
    g1mshort[:, 2:6, :] .= g1m[:, 4:8, :]
    g2mshort[:, 2:6, :] .= g2m[:, 4:8, :]
    p1 = DrugResponseModel.plot_fig1(concs[:, 2], G1short[:, :, 2], g1mshort[:, :, 2, 1], "Dynamical Model Fits - Doxorubicin", "G1", "A", :PuBu_6)
    p2 = DrugResponseModel.plot_fig1(concs[:, 2], G2short[:, :, 2], g2mshort[:, :, 2, 1], "Dynamical Model Fits - Doxorubicin", "S/G2", "B", :PuBu_6)
    p5 = DrugResponseModel.plot_fig1(concs[:, 4], G1short[:, :, 4], g1mshort[:, :, 4, 1], "Dynamical Model Fits - Paclitaxel", "G1", "C", :PuBu_6)
    p6 = DrugResponseModel.plot_fig1(concs[:, 4], G2short[:, :, 4], g2mshort[:, :, 4, 1], "Dynamical Model Fits - Paclitaxel", "S/G2", "D", :PuBu_6)
    p9 = DrugResponseModel.plot_fig1(concs[:, 5], G1short[:, :, 5], g1mshort[:, :, 5, 1], "Dynamical Model Fits - Palbociclib", "G1", "E", :PuBu_6)
    p10 = DrugResponseModel.plot_fig1(concs[:, 5], G2short[:, :, 5], g2mshort[:, :, 5, 1], "Dynamical Model Fits - Palbociclib", "S/G2", "F", :PuBu_6)
    p3 = DrugResponseModel.plot_pG1(duration_efcs[:, :, 3], 60, "Gemcitabine", "Avg phase duration [hr]", "I", 0.35)
    p4 = DrugResponseModel.plot_pG1(death_efcs[:, :, 3], 0.2, "Gemcitabine", "Cell death probability", "J", 0.06)
    p7 = DrugResponseModel.plot_pG1(duration_efcs[:, :, 5], 60, "Paclitaxel", "Avg phase duration [hr]", "K", 0.35)
    p8 = DrugResponseModel.plot_pG1(death_efcs[:, :, 5], 0.2, "Paclitaxel", "Cell death probability", "L", 0.06)
    p11 = DrugResponseModel.plot_pG1(duration_efcs[:, :, 1], 60, "Lapatinib", "Avg phase duration [hr]", "G", 0.35)
    p12 = DrugResponseModel.plot_pG1(death_efcs[:, :, 1], 0.2, "Lapatinib", "Cell death probability", "H", 0.06)
    figureS2 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (1700, 1100), layout = (3, 4))
    savefig(figureS2, "figureS2.svg")
end