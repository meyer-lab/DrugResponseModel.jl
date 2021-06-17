""" Figure 4: drug combination. """


function SSEs_combination()
    cellnum = DrugResponseModel.output_Bliss_cellnum() # bliss on cell number [189, 5, 5, 10]
    g = JLD.load("data/GC.jld")["GC"] # experimental data
    gT = g[1, 1:189, :, :] .+ g[2, 1:189, :, :] # [189, 6, 8]

    p = [
        0.151871,
        0.677644,
        1.24172,
        1.22174,
        0.0796602,
        1.65385,
        0.000103771,
        0.000137389,
        0.00010086,
        0.000107739,
        0.000100047,
        0.0496591,
        0.243334,
        0.61755,
        0.404729,
        0.566801,
        1.04495,
        1.81665,
        0.000103482,
        0.000114637,
        0.000101504,
        0.000185468,
        0.000101091,
        0.199953,
        63.3619,
        2.91267,
        0.000296213,
        0.0350896,
        1.37231,
        0.0599936,
        2.22239,
        0.683019,
        0.00978513,
        0.000100827,
        0.0419537,
        0.000107985,
        0.000116839,
        0.0347002,
        8.62214,
        19.5609,
        3.99883,
        0.187608,
        0.431694,
        0.0710716,
        3.17961,
        0.550232,
        0.00308626,
        0.000101195,
        0.000102833,
        0.0066353,
        0.0241068,
        0.000104833,
        22.9465,
        0.621644,
        3.05299,
        0.156228,
        0.361387,
        1.03003,
        0.229558,
        0.807047,
        0.000492391,
        0.000100816,
        0.000102552,
        0.0165746,
        0.000136098,
        0.0587433,
        3.9392,
        0.228883,
        1.10867,
        0.24115,
        2.52559,
        0.288521,
    ]
    combinEffects = DrugResponseModel.my_helper(p)

    Gsim = zeros(2, 189, 6, 8)
    t = LinRange(0.0, 95.0, 189)

    for j = 1:8
        for i = 1:6 # concentration number
            Gsim[1, :, i, j], Gsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
        end
    end
    GsimT = Gsim[1, :, :, :] .+ Gsim[2, :, :, :] # Total model predictions    [189, 6, 8]

    SSEs = zeros(7, 2) # dim1: Bliss on cell number - exp, dim2: Bliss on model - exp

    # avg(norm(Bliss on cell number - exp))
    SSEs[1, 1] = sum((cellnum[1:189, :, 3, 4] .- gT[1:189, 2:6, 1])) / (5 * 189)
    SSEs[2, 1] = sum((cellnum[1:189, :, 3, 9] .- gT[1:189, 2:6, 2])) / (5 * 189)
    SSEs[3, 1] = sum((cellnum[1:189, 3, :, 9] .- gT[1:189, 2:6, 3])) / (5 * 189)
    SSEs[4, 1] = sum((cellnum[1:189, :, 3, 2] .- gT[1:189, 2:6, 4])) / (5 * 189)
    SSEs[5, 1] = sum((cellnum[1:189, 4, :, 4] .- gT[1:189, 2:6, 5])) / (5 * 189)
    SSEs[6, 1] = sum((cellnum[1:189, 4, :, 2] .- gT[1:189, 2:6, 6])) / (5 * 189)
    SSEs[7, 1] = sum((cellnum[1:189, 2, :, 5] .- gT[1:189, 2:6, 7])) / (5 * 189)

    # avg(norm(model predictions - exp))
    for i = 1:7
        SSEs[i, 2] = sum((GsimT[1:189, 2:6, i] .- gT[1:189, 2:6, i])) / (5 * 189)
    end

    return SSEs
end

function plot_sse_sep(nams, SSE)
    ctg = repeat(["Blisscell number - exp", "Model pred - exp"], inner = 5)
    p = groupedbar(
        SSE,
        group = ctg,
        ylabel = "Cell number difference",
        title = "Comparison",
        bar_width = 0.45,
        xticks = (1:10, nams),
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
    p
end

function plotSSEs()
    SSE = SSEs_combination()
    ctg = repeat(["cell number", "model"], inner = 7)
    nams = repeat(
        [
            "Palb 50 nM + Lpts",
            "Palb 50 nM + Gems",
            "Gem 10 nM + Palbs",
            "Gem 10 nM + Lpts",
            "Lpt 100 nM + Gems",
            "Lpt 100 nM + Palbs",
            "Dox 20 nM + Gems",
        ],
        outer = 2,
    )
    p = groupedbar(
        SSE,
        group = ctg,
        yrotation = 60,
        ylabel = "Combinations",
        xlabel = "SSE",
        yflip = true,
        title = "Sum of Squared Errors",
        bar_width = 0.45,
        yticks = (1:14, nams),
        lw = 0,
        framestyle = :box,
        orientation = :horizontal,
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
    annotate!(-1, 0.0, text("j", :black, :left, Plots.font("Helvetica Bold", 15)))
    p
end

function figure4()

    p = plotSSEs()
    savefig(p, "figure4.svg")
end

function figure4a()
    cellnum = DrugResponseModel.output_Bliss_cellnum() # bliss on cell number [189, 5, 5, 10]
    g = JLD.load("g.jld")["g"] # experimental data
    gT = g[1, 1:189, :, :] .+ g[2, 1:189, :, :] # [189, 6, 8]

    p = [
        0.151871,
        0.677644,
        1.24172,
        1.22174,
        0.0796602,
        1.65385,
        0.000103771,
        0.000137389,
        0.00010086,
        0.000107739,
        0.000100047,
        0.0496591,
        0.243334,
        0.61755,
        0.404729,
        0.566801,
        1.04495,
        1.81665,
        0.000103482,
        0.000114637,
        0.000101504,
        0.000185468,
        0.000101091,
        0.199953,
        63.3619,
        2.91267,
        0.000296213,
        0.0350896,
        1.37231,
        0.0599936,
        2.22239,
        0.683019,
        0.00978513,
        0.000100827,
        0.0419537,
        0.000107985,
        0.000116839,
        0.0347002,
        8.62214,
        19.5609,
        3.99883,
        0.187608,
        0.431694,
        0.0710716,
        3.17961,
        0.550232,
        0.00308626,
        0.000101195,
        0.000102833,
        0.0066353,
        0.0241068,
        0.000104833,
        22.9465,
        0.621644,
        3.05299,
        0.156228,
        0.361387,
        1.03003,
        0.229558,
        0.807047,
        0.000492391,
        0.000100816,
        0.000102552,
        0.0165746,
        0.000136098,
        0.0587433,
        3.9392,
        0.228883,
        1.10867,
        0.24115,
        2.52559,
        0.288521,
    ]
    combinEffects = DrugResponseModel.my_helper(p)

    Gsim = zeros(2, 189, 6, 8)
    t = LinRange(0.0, 95.0, 189)

    for j = 1:8
        for i = 1:6 # concentration number
            Gsim[1, :, i, j], Gsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
        end
    end
    GsimT = Gsim[1, :, :, :] .+ Gsim[2, :, :, :] # Total model predictions    [189, 6, 8]

    SSEs = zeros(5, 2, 7) # dim1: Bliss on cell number - exp, dim2: Bliss on model - exp

    # avg(norm(Bliss on cell number - exp))
    SSEs[:, 1, 1] = sum((cellnum[1:189, :, 3, 4] .- gT[1:189, 2:6, 1]), dims = 1) / 189
    SSEs[:, 1, 2] = sum((cellnum[1:189, :, 3, 9] .- gT[1:189, 2:6, 2]), dims = 1) / 189
    SSEs[:, 1, 3] = sum((cellnum[1:189, 3, :, 9] .- gT[1:189, 2:6, 3]), dims = 1) / 189
    SSEs[:, 1, 4] = sum((cellnum[1:189, :, 3, 2] .- gT[1:189, 2:6, 4]), dims = 1) / 189
    SSEs[:, 1, 5] = sum((cellnum[1:189, 4, :, 4] .- gT[1:189, 2:6, 5]), dims = 1) / 189
    SSEs[:, 1, 6] = sum((cellnum[1:189, 4, :, 2] .- gT[1:189, 2:6, 6]), dims = 1) / 189
    SSEs[:, 1, 7] = sum((cellnum[1:189, 2, :, 5] .- gT[1:189, 2:6, 7]), dims = 1) / 189

    # avg(norm(model predictions - exp))
    for i = 1:7
        SSEs[:, 2, i] = sum((GsimT[1:189, 2:6, i] .- gT[1:189, 2:6, i]), dims = 1) / 189
    end

    p1 = plot_sse_sep(
        repeat(
            ["palbo 50 nM", "palbo 50 nM + lap25 nM", "palbo 50 nM + lap 50 nM", "palbo 50 nM + lap 100 nM", "palbo 50 nM + lap 250 nM"],
            outer = 2,
        ),
        SSEs[:, :, 1],
    )
    p2 = plot_sse_sep(
        repeat(["palbo 50 nM", "palbo 50 nM + gem 5 nM", "palbo 50 nM + gem 10 nM", "palbo 50 nM + gem 17 nM", "palbo 50 nM + gem 30 nM"], outer = 2),
        SSEs[:, :, 2],
    )
    p3 = plot_sse_sep(
        repeat(
            ["gem 10 nM", "gem 10 nM + palbo 25 nM", "gem 10 nM + palbo 50 nM", "gem 10 nM + palbo 100 nM", "gem 10 nM + palbo 250 nM"],
            outer = 2,
        ),
        SSEs[:, :, 3],
    )
    p4 = plot_sse_sep(
        repeat(["gem 10 nM", "gem 10 nM + lap 25 nM", "gem 10 nM + lap 50 nM", "gem 10 nM + lap 100 nM", "gem 10 nM + lap 250 nM"], outer = 2),
        SSEs[:, :, 4],
    )
    p5 = plot_sse_sep(
        repeat(
            ["lap 100 nM", "lap 100 nM + palbo 25 nM", "lap 100 nM + palbo 50 nM", "lap 100 nM + palbo 100 nM", "lap 100 nM + palbo 250 nM"],
            outer = 2,
        ),
        SSEs[:, :, 5],
    )
    p6 = plot_sse_sep(
        repeat(["lap 100 nM", "lap 100 nM + gem 5 nM", "lap 100 nM + gem 10 nM", "lap 100 nM + gem 17 nM", "lap 100 nM + gem 30 nM"], outer = 2),
        SSEs[:, :, 6],
    )
    p7 = plot_sse_sep(
        repeat(["dox 20 nM", "dox 20 nM + gem 5 nM", "dox 20 nM + gem 10 nM", "dox 20 nM + gem 17 nM", "dox 20 nM + gem 30 nM"], outer = 2),
        SSEs[:, :, 7],
    )
    p8 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, size = (2000, 700), layout = (2, 4))
    savefig(p, "fig4a.svg")
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

function output_combination()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    t = LinRange(0.0, 95.0, 189)
    p = parameters()
    efcs = getODEparams(p, concs)
    gem17 = DrugResponseModel.find_gem17(p)
    efcs[:, 8, 3] = efcs[:, 7, 3]
    efcs[:, 7, 3] = gem17

    # Bliss on Model
    LPT_DOX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 2])
    LPT_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3])
    LPT_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 4])
    LPT_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5])
    DOX_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 3])
    DOX_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 4])
    DOX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 5])
    GEM_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 4])
    GEM_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 5])
    TAX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 4], efcs[:, :, 5])


    ########### lapatinibs + Doxorubicins
    LapDox = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapDox[1, :, i, j], LapDox[2, :, i, j], _ = predict(LPT_DOX[:, i, j], LPT_DOX[:, 1, 1], t)
        end
    end
    LapDox[3, :, :, :] .= LapDox[1, :, :, :] .+ LapDox[2, :, :, :]

    ########### lapatinibs + Gemcitabines
    LapGem = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapGem[1, :, i, j], LapGem[2, :, i, j], _ = predict(LPT_GEM[:, i, j], LPT_GEM[:, 1, 1], t)
        end
    end
    LapGem[3, :, :, :] .= LapGem[1, :, :, :] .+ LapGem[2, :, :, :]


    ########### lapatinibs + Taxol
    LapTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapTax[1, :, i, j], LapTax[2, :, i, j], _ = predict(LPT_TAX[:, i, j], LPT_TAX[:, 1, 1], t)
        end
    end
    LapTax[3, :, :, :] .= LapTax[1, :, :, :] .+ LapTax[2, :, :, :]


    ########### lapatinibs + Palbociclibs
    LapPalbo = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapPalbo[1, :, i, j], LapPalbo[2, :, i, j], _ = predict(LPT_PLB[:, i, j], LPT_PLB[:, 1, 1], t)
        end
    end
    LapPalbo[3, :, :, :] .= LapPalbo[1, :, :, :] .+ LapPalbo[2, :, :, :]

    ########### doxorubicins + gemcitabines 
    DoxGem = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxGem[1, :, i, j], DoxGem[2, :, i, j], _ = predict(DOX_GEM[:, i, j], DOX_GEM[:, 1, 1], t)
        end
    end
    DoxGem[3, :, :, :] .= DoxGem[1, :, :, :] .+ DoxGem[2, :, :, :]

    ########### doxorubicins + taxols 
    DoxTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxTax[1, :, i, j], DoxTax[2, :, i, j], _ = predict(DOX_TAX[:, i, j], DOX_TAX[:, 1, 1], t)
        end
    end
    DoxTax[3, :, :, :] .= DoxTax[1, :, :, :] .+ DoxTax[2, :, :, :]
    

    ########### doxorubicins + palbos
    DoxPalbo = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxPalbo[1, :, i, j], DoxPalbo[2, :, i, j], _ = predict(DOX_PLB[:, i, j], DOX_PLB[:, 1, 1], t)
        end
    end
    DoxPalbo[3, :, :, :] .= DoxPalbo[1, :, :, :] .+ DoxPalbo[2, :, :, :]

    ########### gemcitabines + taxols
    GemTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            GemTax[1, :, i, j], GemTax[2, :, i, j], _ = predict(GEM_TAX[:, i, j], GEM_TAX[:, 1, 1], t)
        end
    end
    GemTax[3, :, :, :] .= GemTax[1, :, :, :] .+ GemTax[2, :, :, :]

    ########### gemcitabines + palbos
    GemPalb = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            GemPalb[1, :, i, j], GemPalb[2, :, i, j], _ = predict(GEM_PLB[:, i, j], GEM_PLB[:, 1, 1], t)
        end
    end
    GemPalb[3, :, :, :] .= GemPalb[1, :, :, :] .+ GemPalb[2, :, :, :]

    ########### gemcitabines + palbos
    TaxPlb = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            TaxPlb[1, :, i, j], TaxPlb[2, :, i, j], _ = predict(TAX_PLB[:, i, j], TAX_PLB[:, 1, 1], t)
        end
    end
    TaxPlb[3, :, :, :] .= TaxPlb[1, :, :, :] .+ TaxPlb[2, :, :, :]

    # G1
    df1 = DataFrames.DataFrame(
        palbo50c = palbo_alone[1, :, 2],
        palbo50_lpt25 = palbo50Lap[1, :, 1],
        palbo50_lpt50 = palbo50Lap[1, :, 2],
        palbo50_lpt100 = palbo50Lap[1, :, 3],
        palbo50_lpt250 = palbo50Lap[1, :, 4],
        palbo50 = palbo_alone[1, :, 2],
        palbo50_gem5 = palbo50Gem[1, :, 1],
        palbo50_gem10 = palbo50Gem[1, :, 2],
        palbo50_gem17 = palbo50Gem[1, :, 3],
        palbo50_gem30 = palbo50Gem[1, :, 4],
        gem10c = gemcitabine_alone[1, :, 2],
        gem10lpt25 = Gem10Lap[1, :, 1],
        gem10lpt50 = Gem10Lap[1, :, 2],
        gem10lpt100 = Gem10Lap[1, :, 3],
        gem10lpt250 = Gem10Lap[1, :, 4],
        gem10 = gemcitabine_alone[1, :, 2],
        gem10_palbo25 = Gem10palbo[1, :, 1],
        gem10_palbo50 = Gem10palbo[1, :, 2],
        gem10_palbo100 = Gem10palbo[1, :, 3],
        gem10_palbo250 = Gem10palbo[1, :, 4],
        lpt100c = lapatinib_alone[1, :, 3],
        lpt100palbo25 = Lap100palbo[1, :, 1],
        lpt100palbo50 = Lap100palbo[1, :, 2],
        lpt100palbo100 = Lap100palbo[1, :, 3],
        lpt100palbo250 = Lap100palbo[1, :, 4],
        lpt100 = lapatinib_alone[1, :, 3],
        lpt100_gem5 = lap100gem[1, :, 1],
        lpt100_gem10 = lap100gem[1, :, 2],
        lpt100_gem17 = lap100gem[1, :, 3],
        lpt100_gem30 = lap100gem[1, :, 4],
        dox20 = dox_alone[1, :],
        dox20_gem5 = dox20gem[1, :, 1],
        dox20_gem10 = dox20gem[1, :, 2],
        dox20_gem17 = dox20gem[1, :, 3],
        dox20_gem30 = dox20gem[1, :, 4],
        pax2 = pax2_alone[1, :],
        pax2_lpt25 = Pax2_lap[1, :, 1],
        pax2_lpt50 = Pax2_lap[1, :, 2],
        pax2_lpt100 = Pax2_lap[1, :, 3],
        pax2_lpt250 = Pax2_lap[1, :, 4],
    )

    # G2
    df2 = DataFrames.DataFrame(
        palbo50c = palbo_alone[2, :, 2],
        palbo50_lpt25 = palbo50Lap[2, :, 1],
        palbo50_lpt50 = palbo50Lap[2, :, 2],
        palbo50_lpt100 = palbo50Lap[2, :, 3],
        palbo50_lpt250 = palbo50Lap[2, :, 4],
        palbo50 = palbo_alone[2, :, 2],
        palbo50_gem5 = palbo50Gem[2, :, 1],
        palbo50_gem10 = palbo50Gem[2, :, 2],
        palbo50_gem17 = palbo50Gem[2, :, 3],
        palbo50_gem30 = palbo50Gem[2, :, 4],
        gem10c = gemcitabine_alone[2, :, 2],
        gem10lpt25 = Gem10Lap[2, :, 1],
        gem10lpt50 = Gem10Lap[2, :, 2],
        gem10lpt100 = Gem10Lap[2, :, 3],
        gem10lpt250 = Gem10Lap[2, :, 4],
        gem10 = gemcitabine_alone[2, :, 2],
        gem10_palbo25 = Gem10palbo[2, :, 1],
        gem10_palbo50 = Gem10palbo[2, :, 2],
        gem10_palbo100 = Gem10palbo[2, :, 3],
        gem10_palbo250 = Gem10palbo[2, :, 4],
        lpt100 = lapatinib_alone[2, :, 3],
        lpt100palbo25 = Lap100palbo[2, :, 1],
        lpt100palbo50 = Lap100palbo[2, :, 2],
        lpt100palbo100 = Lap100palbo[2, :, 3],
        lpt100palbo250 = Lap100palbo[2, :, 4],
        lpt100c = lapatinib_alone[2, :, 3],
        lpt100_gem5 = lap100gem[2, :, 1],
        lpt100_gem10 = lap100gem[2, :, 2],
        lpt100_gem17 = lap100gem[2, :, 3],
        lpt100_gem30 = lap100gem[2, :, 4],
        dox20 = dox_alone[2, :],
        dox20_gem5 = dox20gem[2, :, 1],
        dox20_gem10 = dox20gem[2, :, 2],
        dox20_gem17 = dox20gem[2, :, 3],
        dox20_gem30 = dox20gem[2, :, 4],
        pax2 = pax2_alone[2, :],
        pax2_lpt25 = Pax2_lap[2, :, 1],
        pax2_lpt50 = Pax2_lap[2, :, 2],
        pax2_lpt100 = Pax2_lap[2, :, 3],
        pax2_lpt250 = Pax2_lap[2, :, 4],
    )

    # total
    df3 = DataFrames.DataFrame(
        palbo50c = palbo_alone[3, :, 2],
        palbo50_lpt25 = palbo50Lap[3, :, 1],
        palbo50_lpt50 = palbo50Lap[3, :, 2],
        palbo50_lpt100 = palbo50Lap[3, :, 3],
        palbo50_lpt250 = palbo50Lap[3, :, 4],
        palbo50 = palbo_alone[3, :, 2],
        palbo50_gem5 = palbo50Gem[3, :, 1],
        palbo50_gem10 = palbo50Gem[3, :, 2],
        palbo50_gem17 = palbo50Gem[3, :, 3],
        palbo50_gem30 = palbo50Gem[3, :, 4],
        gem10c = gemcitabine_alone[3, :, 2],
        gem10lpt25 = Gem10Lap[3, :, 1],
        gem10lpt50 = Gem10Lap[3, :, 2],
        gem10lpt100 = Gem10Lap[3, :, 3],
        gem10lpt250 = Gem10Lap[3, :, 4],
        gem10 = gemcitabine_alone[3, :, 2],
        gem10_palbo25 = Gem10palbo[3, :, 1],
        gem10_palbo50 = Gem10palbo[3, :, 2],
        gem10_palbo100 = Gem10palbo[3, :, 3],
        gem10_palbo250 = Gem10palbo[3, :, 4],
        lpt100c = lapatinib_alone[3, :, 3],
        lpt100palbo25 = Lap100palbo[3, :, 1],
        lpt100palbo50 = Lap100palbo[3, :, 2],
        lpt100palbo100 = Lap100palbo[3, :, 3],
        lpt100palbo250 = Lap100palbo[3, :, 4],
        lpt100 = lapatinib_alone[3, :, 3],
        lpt100_gem5 = lap100gem[3, :, 1],
        lpt100_gem10 = lap100gem[3, :, 2],
        lpt100_gem17 = lap100gem[3, :, 3],
        lpt100_gem30 = lap100gem[3, :, 4],
        dox20 = dox_alone[3, :],
        dox20_gem5 = dox20gem[3, :, 1],
        dox20_gem10 = dox20gem[3, :, 2],
        dox20_gem17 = dox20gem[3, :, 3],
        dox20_gem30 = dox20gem[3, :, 4],
        pax2 = pax2_alone[2, :],
        pax2_lpt25 = Pax2_lap[3, :, 1],
        pax2_lpt50 = Pax2_lap[3, :, 2],
        pax2_lpt100 = Pax2_lap[3, :, 3],
        pax2_lpt250 = Pax2_lap[3, :, 4],
    )

    XLSX.writetable("G1.xlsx", df1)
    XLSX.writetable("G2.xlsx", df2)
    XLSX.writetable("Total.xlsx", df3)
end
