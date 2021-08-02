""" Plot the synergy/antagonism summary plot"""


function plot_summary()
    cn = LinRange(0.0, 250.0, 20)
    conc = hcat(cn, cn, cn)
    p = [
        70.0,
        3.1,
        0.4,
        0.4,
        0.4,
        0.4,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        70.0,
        3.1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0.08,
        0.08,
        0.08,
        0.08,
        0,
        0,
        0,
        0,
        70.0,
        3.1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0.02,
        0.02,
        0.02,
        0.02,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
    ]
    ef = getODEparams(p, conc)

    Pcombin = zeros(12, 20, 2)
    for i = 1:20
        Pcombin[:, i, 1] = DrugResponseModel.Bliss_params_unit(ef[:, 4, 1], ef[:, i, 2], hcat(ef[:, 1, 1], ef[:, 1, 1])) # A & B
        Pcombin[:, i, 2] = DrugResponseModel.Bliss_params_unit(ef[:, 4, 1], ef[:, i, 3], hcat(ef[:, 1, 1], ef[:, 1, 1])) # A & C
    end

    expected_cellnum = zeros(3, 20, 2) # from combination
    for i = 1:20
        expected_cellnum[1, i, 1], expected_cellnum[2, i, 1], _ = DrugResponseModel.predict(Pcombin[:, i, 1], ef[:, 1, 1], 96.0)
        expected_cellnum[1, i, 2], expected_cellnum[2, i, 2], _ = DrugResponseModel.predict(Pcombin[:, i, 2], ef[:, 1, 1], 96.0)
    end
    expected_cellnum[3, :, :] .= expected_cellnum[1, :, :] .+ expected_cellnum[2, :, :]

    # single treatment expected cell numbers
    single_cellnums = zeros(3, 20)
    for i = 1:20
        single_cellnums[1, i], single_cellnums[2, i], _ = predict(ef[:, i, 1], ef[:, 1, 1], 96.0)
    end
    single_cellnums[3, :] .= single_cellnums[1, :] .+ single_cellnums[2, :] # A

    function unit_plt(cnc, ycombin, ysingle, subPlabel)
        plot(
            cnc,
            ycombin ./ ysingle[1, 1],
            label = "A+B",
            xlabel = "concentration [nM]",
            ylabel = "cell #",
            fg_legend = :transparent,
            lw = 3,
            alpha = 0.8,
            titlefont = Plots.font("Helvetica", 14),
            legendfont = Plots.font("Helvetica", 11),
            guidefont = Plots.font("Helvetica", 14),
            xtickfont = Plots.font("Helvetica", 14),
            ytickfont = Plots.font("Helvetica", 14),
            bottom_margin = 1.25cm,
            top_margin = 1.25cm,
            left_margin = 1.25cm,
            right_margin = 1.25cm,
        )
        plot!(cnc, ysingle ./ ysingle[1, 1], alpha = 0.8, lw = 3, label = "A")
        annotate!(-0.1, 2.7, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
        ylims!((-0.05, 1.2))
    end
    p1 = unit_plt(cn, expected_cellnum[3, :, 1], single_cellnums[3, :], "A") # A
    p2 = unit_plt(cn, expected_cellnum[3, :, 2], single_cellnums[3, :], "B") # B
    P = plot(p1, p2, layout = (1, 2), size = (900, 350))
    savefig(P, "summary.svg")
end



""" Plot the synergy/antagonism summary plot"""


function fake_drugs()
    conc = hcat([0.0, 25.0, 50.0, 100.0, 250.0], [0.0, 25.0, 50.0, 100.0, 250.0], [0.0, 25.0, 50.0, 100.0, 250.0])
    p = [
        50.0,
        1.1,
        0.1,
        0.1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        50.0,
        1.1,
        1,
        1,
        1,
        1,
        1,
        1,
        0.02,
        0.02,
        0,
        0,
        0,
        0,
        50.0,
        1.1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0.01,
        0.01,
        0.01,
        0.01,
        1,
        1,
        1,
        1,
        1,
        1,
    ]
    ef = getODEparams(p, conc)

    # single treatment expected cell numbers
    single_cellnums = zeros(3, 5, 3)
    for i = 1:5
        single_cellnums[1, i, 1], single_cellnums[2, i, 1], _ = predict(ef[:, i, 1], ef[:, 1, 1], 96.0)
        single_cellnums[1, i, 2], single_cellnums[2, i, 2], _ = predict(ef[:, i, 2], ef[:, 1, 2], 96.0)
        single_cellnums[1, i, 3], single_cellnums[2, i, 3], _ = predict(ef[:, i, 3], ef[:, 1, 3], 96.0)
    end
    single_cellnums[3, :, :] .= single_cellnums[1, :, :] .+ single_cellnums[2, :, :] # total

    # combination treatment expected cell numbers
    Pcombin = zeros(10, 12, 5) # first 5 drug A + drug B, last 5: drugA + drugC
    for j = 1:5
        for i = 1:5
            Pcombin[j, :, i] = DrugResponseModel.Bliss_params_unit(ef[:, j, 1], ef[:, i, 2], hcat(ef[:, 1, 1], ef[:, 1, 1]))
            Pcombin[j + 5, :, i] = DrugResponseModel.Bliss_params_unit(ef[:, j, 1], ef[:, i, 3], hcat(ef[:, 1, 1], ef[:, 1, 1]))
        end
    end

    expected_cellnums = zeros(3, 10, 5)

    for j = 1:5
        for i = 1:5
            expected_cellnums[1, j, i], expected_cellnums[2, j, i], _ = predict(Pcombin[j, :, i], ef[:, 1, 1], 96.0)
            expected_cellnums[1, j + 5, i], expected_cellnums[2, j + 5, i], _ = predict(Pcombin[j + 5, :, i], ef[:, 1, 1], 96.0)
        end
    end
    expected_cellnums[3, :, :] .= expected_cellnums[1, :, :] .+ expected_cellnums[2, :, :] # total

    onee = expected_cellnums[3, 1:5, :]
    twoo = expected_cellnums[3, 6:10, :]
    function unit_plt(cnc, onee, subPlabel, d, titl)
        plot(
            cnc,
            onee[:, 1],
            label = "drugA",
            ylabel = "cell #",
            fg_legend = :transparent,
            xlabel = "concentration [nM]",
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
            title = titl,
        )
        plot!(cnc, onee[:, 2], label = "drugA+drug$d 25", lw = 3, alpha = 0.8)
        plot!(cnc, onee[:, 3], label = "drugA+drug$d 50", lw = 3, alpha = 0.8)
        plot!(cnc, onee[:, 4], label = "drugA+drug$d 100", lw = 3, alpha = 0.8)
        plot!(cnc, onee[:, 5], label = "drugA+drug $d 250", lw = 3, alpha = 0.8)
        # annotate!(-0.1, 2.7, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
        ylims!((-0.05, 11.0))
    end

    p1 = unit_plt(conc[:, 1], onee, "A", "B", "G1 cell arrest and death")
    p2 = unit_plt(conc[:, 1], twoo, "B", "C", "G1 cell arrest S/G2 death")
    p = plot(p1, p2, layout = (1, 2), size = (900, 400))
    savefig(p, "summ.svg")
end
