""" Plot the synergy/antagonism summary plot"""

function find_gem17(p)
    # Interpolation to find the parameters for 17 nM.
    hill(p, c) = p[2] + (p[3] - p[2]) / (1 + ((p[1] / c)^p[4]))
    gemc_hillParams = zeros(12, 4) # [a1, a2, b1, b2, b3, b4, d1, d2, d3, d4, d5, d6] x [EC50, min, max, k]
    gemc_hillParams[:, 1] .= p[15] # ec50
    gemc_hillParams[:, 4] .= p[16] # k
    gemc_hillParams[1:6, 2] = p[71:76]
    gemc_hillParams[7:12, 2] .= 0.0
    gemc_hillParams[:, 3] .= p[17:28]
    GEM17 = zeros(12)
    for i = 1:length(GEM17)
        GEM17[i] = hill(gemc_hillParams[i, :], 17.0)
    end
    GEM17
end

function plot_summary()
    # load experimental data
    g = JLD.load("data/GC.jld")["GC"]

    concs, _, g1s1, g2s1 = load(189, 1)
    _, _, g1s2, g2s2 = load(189, 2)
    _, _, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    gm = mean(g1S, dims = 4) .+ mean(g2S, dims = 4)

    # load bliss on cell numbers
    bliss = DrugResponseModel.blissCellNum(gm[189, :, :, 1]; n = 8)

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
    efcs[:, 8, 3] = efcs[:, 7, 3]
    efcs[:, 7, 3] = find_gem17(ps)

    # combination effects:
    lpt_palbo = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5]; n = 8) # A: palbo 50 + lapatinib [25, 50, 100, 250]
    gem_palbo = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 5]; n = 8) # B: palbo 50 + gemcitabine [5, 10, 17, 30]
    dox_gem = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 3]; n = 8) # C: dox 20 + gemcitabine [5, 10, 17, 30]
    lap_gem = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3]; n = 8) # D: lap 100 + gemcitabine [5, 10, 17, 30]

    expected_cellnum = zeros(3, 5, 4)
    for i = 4:7
        expected_cellnum[1, i - 2, 1], expected_cellnum[2, i - 2, 1], _ = DrugResponseModel.predict(lpt_palbo[:, i, 5], efcs[:, 1, 1], 96.0)
        expected_cellnum[1, i - 2, 2], expected_cellnum[2, i - 2, 2], _ = DrugResponseModel.predict(gem_palbo[:, i + 1, 5], efcs[:, 1, 1], 96.0)
        expected_cellnum[1, i - 2, 3], expected_cellnum[2, i - 2, 3], _ = DrugResponseModel.predict(dox_gem[:, 4, i + 1], efcs[:, 1, 1], 96.0)
        expected_cellnum[1, i - 2, 4], expected_cellnum[2, i - 2, 4], _ = DrugResponseModel.predict(lap_gem[:, 6, i], efcs[:, 1, 1], 96.0)
    end
    expected_cellnum[1, 1, 1], expected_cellnum[2, 1, 1], _ = DrugResponseModel.predict(lpt_palbo[:, 1, 5], efcs[:, 1, 1], 96.0)
    expected_cellnum[1, 1, 2], expected_cellnum[2, 1, 2], _ = DrugResponseModel.predict(gem_palbo[:, 1, 5], efcs[:, 1, 1], 96.0)
    expected_cellnum[1, 1, 3], expected_cellnum[2, 1, 3], _ = DrugResponseModel.predict(dox_gem[:, 4, 1], efcs[:, 1, 1], 96.0)
    expected_cellnum[1, 1, 4], expected_cellnum[2, 1, 4], _ = DrugResponseModel.predict(lap_gem[:, 6, 1], efcs[:, 1, 1], 96.0)

    expected_cellnum[3, :, :] .= expected_cellnum[1, :, :] .+ expected_cellnum[2, :, :]

    function unit_plt(cnc, yexpect, yexp, ybase, subPlabel)
        plot(
            cnc,
            yexpect,
            label = "Model",
            xlabel = "concentration [nM]",
            ylabel = "cell #",
            fg_legend = :transparent,
            lw = 3,
            alpha = 0.8,
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
        plot!(cnc, yexp, alpha = 0.8, lw = 3, label = "Experiment")
        plot!(cnc, ybase, alpha = 0.8, lw = 3, label = "Bliss")
        annotate!(-0.1, 2.7, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
        ylims!((-0.05, 2.5))
    end
    p1 = unit_plt(concs[[1, 4, 5, 6, 7], 1], expected_cellnum[3, :, 1], g[3, 193, 2:6, 1], bliss[[1, 4, 5, 6, 7], 5, 4], "A") # A
    p2 = unit_plt([0, 5, 10, 17, 30], expected_cellnum[3, :, 2], g[3, 193, 2:6, 2], bliss[[1, 5, 6, 7, 8], 5, 9], "B") # B
    p3 = unit_plt([0, 5, 10, 17, 30], expected_cellnum[3, :, 3], g[3, 193, 2:6, 7], bliss[4, [1, 5, 6, 7, 8], 5], "C") # C
    p4 = unit_plt([0, 5, 10, 17, 30], expected_cellnum[3, :, 4], g[3, 193, 2:6, 6], bliss[6, [1, 5, 6, 7, 8], 2], "D") # D
    P = plot(p1, p2, p3, p4, layout = (2, 2), size = (900, 700))
    savefig(P, "summary.svg")
end
