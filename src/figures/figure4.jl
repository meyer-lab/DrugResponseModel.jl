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

    ########### taxol + palbos
    TaxPlb = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            TaxPlb[1, :, i, j], TaxPlb[2, :, i, j], _ = predict(TAX_PLB[:, i, j], TAX_PLB[:, 1, 1], t)
        end
    end
    TaxPlb[3, :, :, :] .= TaxPlb[1, :, :, :] .+ TaxPlb[2, :, :, :]

    # total
    df1 = DataFrames.DataFrame(
        control_control = LPT_DOX[3, :, 1, 1],
        dox1 = LPT_DOX[3, :, 1, 2],
        dox10 = LPT_DOX[3, :, 1, 3],
        dox20 = LPT_DOX[3, :, 1, 4],
        dox50 = LPT_DOX[3, :, 1, 5],
        dox100 = LPT_DOX[3, :, 1, 6],
        dox250 = LPT_DOX[3, :, 1, 7],
        dox500 = LPT_DOX[3, :, 1, 8],
        lpt5 = LPT_DOX[3, :, 2, 1],
        lpt5dox1 = LPT_DOX[3, :, 2, 2],
        lpt5dox10 = LPT_DOX[3, :, 2, 3],
        lpt5dox20 = LPT_DOX[3, :, 2, 4],
        lpt5dox50 = LPT_DOX[3, :, 2, 5],
        lpt5dox100 = LPT_DOX[3, :, 2, 6],
        lpt5dox250 = LPT_DOX[3, :, 2, 7],
        lpt5dox500 = LPT_DOX[3, :, 2, 8],
        lpt10 = LPT_DOX[3, :, 3, 1],
        lpt10dox1 = LPT_DOX[3, :, 3, 2],
        lpt10dox10 = LPT_DOX[3, :, 3, 3],
        lpt10dox20 = LPT_DOX[3, :, 3, 4],
        lpt10dox50 = LPT_DOX[3, :, 3, 5],
        lpt10dox100 = LPT_DOX[3, :, 3, 6],
        lpt10dox250 = LPT_DOX[3, :, 3, 7],
        lpt10dox500 = LPT_DOX[3, :, 3, 8],
        lpt25 = LPT_DOX[3, :, 4, 1],
        lpt25dox1 = LPT_DOX[3, :, 4, 2],
        lpt25dox10 = LPT_DOX[3, :, 4, 3],
        lpt25dox20 = LPT_DOX[3, :, 4, 4],
        lpt25dox50 = LPT_DOX[3, :, 4, 5],
        lpt25dox100 = LPT_DOX[3, :, 4, 6],
        lpt25dox250 = LPT_DOX[3, :, 4, 7],
        lpt25dox500 = LPT_DOX[3, :, 4, 8],
        lpt50 = LPT_DOX[3, :, 5, 1],
        lpt50dox1 = LPT_DOX[3, :, 5, 2],
        lpt50dox10 = LPT_DOX[3, :, 5, 3],
        lpt50dox20 = LPT_DOX[3, :, 5, 4],
        lpt50dox50 = LPT_DOX[3, :, 5, 5],
        lpt50dox100 = LPT_DOX[3, :, 5, 6],
        lpt50dox250 = LPT_DOX[3, :, 5, 7],
        lpt50dox500 = LPT_DOX[3, :, 5, 8],
        lpt100 = LPT_DOX[3, :, 6, 1],
        lpt100dox1 = LPT_DOX[3, :, 6, 2],
        lpt100dox10 = LPT_DOX[3, :, 6, 3],
        lpt100dox20 = LPT_DOX[3, :, 6, 4],
        lpt100dox50 = LPT_DOX[3, :, 6, 5],
        lpt100dox100 = LPT_DOX[3, :, 6, 6],
        lpt100dox250 = LPT_DOX[3, :, 6, 7],
        lpt100dox500 = LPT_DOX[3, :, 6, 8],
        lpt250 = LPT_DOX[3, :, 7, 1],
        lpt250dox1 = LPT_DOX[3, :, 7, 2],
        lpt250dox10 = LPT_DOX[3, :, 7, 3],
        lpt250dox20 = LPT_DOX[3, :, 7, 4],
        lpt250dox50 = LPT_DOX[3, :, 7, 5],
        lpt250dox100 = LPT_DOX[3, :, 7, 6],
        lpt250dox250 = LPT_DOX[3, :, 7, 7],
        lpt250dox500 = LPT_DOX[3, :, 7, 8],
        lpt500 = LPT_DOX[3, :, 8, 1],
        lpt500dox1 = LPT_DOX[3, :, 8, 2],
        lpt500dox10 = LPT_DOX[3, :, 8, 3],
        lpt500dox20 = LPT_DOX[3, :, 8, 4],
        lpt500dox50 = LPT_DOX[3, :, 8, 5],
        lpt500dox100 = LPT_DOX[3, :, 8, 6],
        lpt500dox250 = LPT_DOX[3, :, 8, 7],
        lpt500dox500 = LPT_DOX[3, :, 8, 8],
    )

    XLSX.writetable("LPT_DOX.xlsx", df1)

    df2 = DataFrames.DataFrame(
        control_control = LPT_GEM[3, :, 1, 1],
        gem0_25 = LPT_GEM[3, :, 1, 2],
        gem1 = LPT_GEM[3, :, 1, 3],
        gem2_5 = LPT_GEM[3, :, 1, 4],
        gem5 = LPT_GEM[3, :, 1, 5],
        gem10 = LPT_GEM[3, :, 1, 6],
        gem17 = LPT_GEM[3, :, 1, 7],
        gem30 = LPT_GEM[3, :, 1, 8],
        lpt5gem0_25 = LPT_GEM[3, :, 2, 2],
        lpt5gem1 = LPT_GEM[3, :, 2, 3],
        lpt5gem2_5 = LPT_GEM[3, :, 2, 4],
        lpt5gem5 = LPT_GEM[3, :, 2, 5],
        lpt5gem10 = LPT_GEM[3, :, 2, 6],
        lpt5gem17 = LPT_GEM[3, :, 2, 7],
        lpt5gem30 = LPT_GEM[3, :, 2, 8],
        lpt10gem0_25 = LPT_GEM[3, :, 3, 2],
        lpt10gem1 = LPT_GEM[3, :, 3, 3],
        lpt10gem2_5 = LPT_GEM[3, :, 3, 4],
        lpt10gem5 = LPT_GEM[3, :, 3, 5],
        lpt10gem10 = LPT_GEM[3, :, 3, 6],
        lpt10gem17 = LPT_GEM[3, :, 3, 7],
        lpt10gem30 = LPT_GEM[3, :, 3, 8],
        lpt25gem0_25 = LPT_GEM[3, :, 4, 2],
        lpt25gem1 = LPT_GEM[3, :, 4, 3],
        lpt25gem2_5 = LPT_GEM[3, :, 4, 4],
        lpt25gem5 = LPT_GEM[3, :, 4, 5],
        lpt25gem10 = LPT_GEM[3, :, 4, 6],
        lpt25gem17 = LPT_GEM[3, :, 4, 7],
        lpt25gem30 = LPT_GEM[3, :, 4, 8],
        lpt50gem0_25 = LPT_GEM[3, :, 5, 2],
        lpt50gem1 = LPT_GEM[3, :, 5, 3],
        lpt50gem2_5 = LPT_GEM[3, :, 5, 4],
        lpt50gem5 = LPT_GEM[3, :, 5, 5],
        lpt50gem10 = LPT_GEM[3, :, 5, 6],
        lpt50gem17 = LPT_GEM[3, :, 5, 7],
        lpt50gem30 = LPT_GEM[3, :, 5, 8],
        lpt100gem0_25 = LPT_GEM[3, :, 6, 2],
        lpt100gem1 = LPT_GEM[3, :, 6, 3],
        lpt100gem2_5 = LPT_GEM[3, :, 6, 4],
        lpt100gem5 = LPT_GEM[3, :, 6, 5],
        lpt100gem10 = LPT_GEM[3, :, 6, 6],
        lpt100gem17 = LPT_GEM[3, :, 6, 7],
        lpt100gem30 = LPT_GEM[3, :, 6, 8],
        lpt250gem0_25 = LPT_GEM[3, :, 7, 2],
        lpt250gem1 = LPT_GEM[3, :, 7, 3],
        lpt250gem2_5 = LPT_GEM[3, :, 7, 4],
        lpt250gem5 = LPT_GEM[3, :, 7, 5],
        lpt250gem10 = LPT_GEM[3, :, 7, 6],
        lpt250gem17 = LPT_GEM[3, :, 7, 7],
        lpt250gem30 = LPT_GEM[3, :, 7, 8],
        lpt500gem0_25 = LPT_GEM[3, :, 8, 2],
        lpt500gem1 = LPT_GEM[3, :, 8, 3],
        lpt500gem2_5 = LPT_GEM[3, :, 8, 4],
        lpt500gem5 = LPT_GEM[3, :, 8, 5],
        lpt500gem10 = LPT_GEM[3, :, 8, 6],
        lpt500gem17 = LPT_GEM[3, :, 8, 7],
        lpt500gem30 = LPT_GEM[3, :, 8, 8],
    )

    XLSX.writetable("LPT_GEM.xlsx", df2)

    df3 = DataFrames.DataFrame(
        control_control = LPT_TAX[3, :, 1, 1],
        taxol0_1 = LPT_TAX[3, :, 1, 2],
        taxol1 = LPT_TAX[3, :, 1, 3],
        taxol2 = LPT_TAX[3, :, 1, 4],
        taxol3 = LPT_TAX[3, :, 1, 5],
        taxol5 = LPT_TAX[3, :, 1, 6],
        taxol7_5 = LPT_TAX[3, :, 1, 7],
        taxol15 = LPT_TAX[3, :, 1, 8],
        lpt5taxol0_1 = LPT_TAX[3, :, 2, 2],
        lpt5taxol1 = LPT_TAX[3, :, 2, 3],
        lpt5taxol2 = LPT_TAX[3, :, 2, 4],
        lpt5taxol3 = LPT_TAX[3, :, 2, 5],
        lpt5taxol5 = LPT_TAX[3, :, 2, 6],
        lpt5taxol7_5 = LPT_TAX[3, :, 2, 7],
        lpt5taxol15 = LPT_TAX[3, :, 2, 8],
        lpt10taxol0_1 = LPT_TAX[3, :, 3, 2],
        lpt10taxol1 = LPT_TAX[3, :, 3, 3],
        lpt10taxol2 = LPT_TAX[3, :, 3, 4],
        lpt10taxol3 = LPT_TAX[3, :, 3, 5],
        lpt10taxol5 = LPT_TAX[3, :, 3, 6],
        lpt10taxol7_5 = LPT_TAX[3, :, 3, 7],
        lpt10taxol15 = LPT_TAX[3, :, 3, 8],
        lpt25taxol0_1 = LPT_TAX[3, :, 4, 2],
        lpt25taxol1 = LPT_TAX[3, :, 4, 3],
        lpt25taxol2 = LPT_TAX[3, :, 4, 4],
        lpt25taxol3 = LPT_TAX[3, :, 4, 5],
        lpt25taxol5 = LPT_TAX[3, :, 4, 6],
        lpt25taxol7_5 = LPT_TAX[3, :, 4, 7],
        lpt25taxol15 = LPT_TAX[3, :, 4, 8],
        lpt50taxol0_1 = LPT_TAX[3, :, 5, 2],
        lpt50taxol1 = LPT_TAX[3, :, 5, 3],
        lpt50taxol2 = LPT_TAX[3, :, 5, 4],
        lpt50taxol3 = LPT_TAX[3, :, 5, 5],
        lpt50taxol5 = LPT_TAX[3, :, 5, 6],
        lpt50taxol7_5 = LPT_TAX[3, :, 5, 7],
        lpt50taxol15 = LPT_TAX[3, :, 5, 8],
        lpt100taxol0_1 = LPT_TAX[3, :, 6, 2],
        lpt100taxol1 = LPT_TAX[3, :, 6, 3],
        lpt100taxol2 = LPT_TAX[3, :, 6, 4],
        lpt100taxol3 = LPT_TAX[3, :, 6, 5],
        lpt100taxol5 = LPT_TAX[3, :, 6, 6],
        lpt100taxol7_5 = LPT_TAX[3, :, 6, 7],
        lpt100taxol15 = LPT_TAX[3, :, 6, 8],
        lpt250taxol0_1 = LPT_TAX[3, :, 7, 2],
        lpt250taxol1 = LPT_TAX[3, :, 7, 3],
        lpt250taxol2 = LPT_TAX[3, :, 7, 4],
        lpt250taxol3 = LPT_TAX[3, :, 7, 5],
        lpt250taxol5 = LPT_TAX[3, :, 7, 6],
        lpt250taxol7_5 = LPT_TAX[3, :, 7, 7],
        lpt250taxol15 = LPT_TAX[3, :, 7, 8],
        lpt500taxol0_1 = LPT_TAX[3, :, 8, 2],
        lpt500taxol1 = LPT_TAX[3, :, 8, 3],
        lpt500taxol2 = LPT_TAX[3, :, 8, 4],
        lpt500taxol3 = LPT_TAX[3, :, 8, 5],
        lpt500taxol5 = LPT_TAX[3, :, 8, 6],
        lpt500taxol7_5 = LPT_TAX[3, :, 8, 7],
        lpt500taxol15 = LPT_TAX[3, :, 8, 8],
    )

    XLSX.writetable("LPT_TAX.xlsx", df3)

    df4 = DataFrames.DataFrame(
        control_control = LPT_PLB[3, :, 1, 1],
        palbo5 = LPT_PLB[3, :, 1, 2],
        palbo10 = LPT_PLB[3, :, 1, 3],
        palbo25 = LPT_PLB[3, :, 1, 4],
        palbo50 = LPT_PLB[3, :, 1, 5],
        palbo100 = LPT_PLB[3, :, 1, 6],
        palbo250 = LPT_PLB[3, :, 1, 7],
        palbo500 = LPT_PLB[3, :, 1, 8],
        lpt5palbo5 = LPT_PLB[3, :, 2, 2],
        lpt5palbo10 = LPT_PLB[3, :, 2, 3],
        lpt5palbo25 = LPT_PLB[3, :, 2, 4],
        lpt5palbo50 = LPT_PLB[3, :, 2, 5],
        lpt5palbo100 = LPT_PLB[3, :, 2, 6],
        lptpalbo250 = LPT_PLB[3, :, 2, 7],
        lpt5palbo500 = LPT_PLB[3, :, 2, 8],
        lpt10palbo5 = LPT_PLB[3, :, 3, 2],
        lpt10palbo10 = LPT_PLB[3, :, 3, 3],
        lpt10palbo25 = LPT_PLB[3, :, 3, 4],
        lpt10palbo50 = LPT_PLB[3, :, 3, 5],
        lpt10palbo100 = LPT_PLB[3, :, 3, 6],
        lpt10palbo250 = LPT_PLB[3, :, 3, 7],
        lpt10palbo500 = LPT_PLB[3, :, 3, 8],
        lpt25palbo5 = LPT_PLB[3, :, 4, 2],
        lpt25palbo10 = LPT_PLB[3, :, 4, 3],
        lpt25palbo25 = LPT_PLB[3, :, 4, 4],
        lpt25palbo50 = LPT_PLB[3, :, 4, 5],
        lpt25palbo100 = LPT_PLB[3, :, 4, 6],
        lpt25palbo250 = LPT_PLB[3, :, 4, 7],
        lpt25palbo500 = LPT_PLB[3, :, 4, 8],
        lpt50palbo5 = LPT_PLB[3, :, 5, 2],
        lpt50palbo10 = LPT_PLB[3, :, 5, 3],
        lpt50palbo25 = LPT_PLB[3, :, 5, 4],
        lpt50palbo50 = LPT_PLB[3, :, 5, 5],
        lpt50palbo100 = LPT_PLB[3, :, 5, 6],
        lpt50palbo250 = LPT_PLB[3, :, 5, 7],
        lpt50palbo500 = LPT_PLB[3, :, 5, 8],
        lpt100palbo5 = LPT_PLB[3, :, 6, 2],
        lpt100palbo10 = LPT_PLB[3, :, 6, 3],
        lpt100palbo25 = LPT_PLB[3, :, 6, 4],
        lpt100palbo50 = LPT_PLB[3, :, 6, 5],
        lpt100palbo100 = LPT_PLB[3, :, 6, 6],
        lpt100palbo250 = LPT_PLB[3, :, 6, 7],
        lpt100palbo500 = LPT_PLB[3, :, 6, 8],
        lpt250palbo5 = LPT_PLB[3, :, 7, 2],
        lpt250palbo10 = LPT_PLB[3, :, 7, 3],
        lpt250palbo25 = LPT_PLB[3, :, 7, 4],
        lpt250palbo50 = LPT_PLB[3, :, 7, 5],
        lpt250palbo100 = LPT_PLB[3, :, 7, 6],
        lpt250palbo250 = LPT_PLB[3, :, 7, 7],
        lpt250palbo500 = LPT_PLB[3, :, 7, 8],
        lpt500palbo5 = LPT_PLB[3, :, 8, 2],
        lpt500palbo10 = LPT_PLB[3, :, 8, 3],
        lpt500palbo25 = LPT_PLB[3, :, 8, 4],
        lpt500palbo50 = LPT_PLB[3, :, 8, 5],
        lpt500palbo100 = LPT_PLB[3, :, 8, 6],
        lpt500palbo250 = LPT_PLB[3, :, 8, 7],
        lpt500palbo500 = LPT_PLB[3, :, 8, 8],
    )

    XLSX.writetable("LPT_PLB.xlsx", df4)

    df5 = DataFrames.DataFrame(
        control_control = DOX_GEM[3, :, 1, 1],
        dox1gem0_25 = DOX_GEM[3, :, 2, 2],
        dox1gem1 = DOX_GEM[3, :, 2, 3],
        dox1gem2_5 = DOX_GEM[3, :, 2, 4],
        dox1gem5 = DOX_GEM[3, :, 2, 5],
        dox1gem10 = DOX_GEM[3, :, 2, 6],
        dox1gem17 = DOX_GEM[3, :, 2, 7],
        dox1gem30 = DOX_GEM[3, :, 2, 8],
        dox10gem0_25 = DOX_GEM[3, :, 3, 2],
        dox10gem1 = DOX_GEM[3, :, 3, 3],
        dox10gem2_5 = DOX_GEM[3, :, 3, 4],
        dox10gem5 = DOX_GEM[3, :, 3, 5],
        dox10gem10 = DOX_GEM[3, :, 3, 6],
        dox10gem17 = DOX_GEM[3, :, 3, 7],
        dox10gem30 = DOX_GEM[3, :, 3, 8],
        dox20gem0_25 = DOX_GEM[3, :, 4, 2],
        dox20gem1 = DOX_GEM[3, :, 4, 3],
        dox20gem2_5 = DOX_GEM[3, :, 4, 4],
        dox20gem5 = DOX_GEM[3, :, 4, 5],
        dox20gem10 = DOX_GEM[3, :, 4, 6],
        dox20gem17 = DOX_GEM[3, :, 4, 7],
        dox20gem30 = DOX_GEM[3, :, 4, 8],
        dox50gem0_25 = DOX_GEM[3, :, 5, 2],
        dox50gem1 = DOX_GEM[3, :, 5, 3],
        dox50gem2_5 = DOX_GEM[3, :, 5, 4],
        dox50gem5 = DOX_GEM[3, :, 5, 5],
        dox50gem10 = DOX_GEM[3, :, 5, 6],
        dox50gem17 = DOX_GEM[3, :, 5, 7],
        dox50gem30 = DOX_GEM[3, :, 5, 8],
        dox125em0_25 = DOX_GEM[3, :, 6, 2],
        dox125gem1 = DOX_GEM[3, :, 6, 3],
        dox125gem2_5 = DOX_GEM[3, :, 6, 4],
        dox125gem5 = DOX_GEM[3, :, 6, 5],
        dox125gem10 = DOX_GEM[3, :, 6, 6],
        dox125gem17 = DOX_GEM[3, :, 6, 7],
        dox125gem30 = DOX_GEM[3, :, 6, 8],
        dox250gem0_25 = DOX_GEM[3, :, 7, 2],
        dox250gem1 = DOX_GEM[3, :, 7, 3],
        dox250gem2_5 = DOX_GEM[3, :, 7, 4],
        dox250gem5 = DOX_GEM[3, :, 7, 5],
        dox250gem10 = DOX_GEM[3, :, 7, 6],
        dox250gem17 = DOX_GEM[3, :, 7, 7],
        dox250gem30 = DOX_GEM[3, :, 7, 8],
        dox500gem0_25 = DOX_GEM[3, :, 8, 2],
        dox500gem1 = DOX_GEM[3, :, 8, 3],
        dox500gem2_5 = DOX_GEM[3, :, 8, 4],
        dox500gem5 = DOX_GEM[3, :, 8, 5],
        dox500gem10 = DOX_GEM[3, :, 8, 6],
        dox500gem17 = DOX_GEM[3, :, 8, 7],
        dox500gem30 = DOX_GEM[3, :, 8, 8],
    )

    XLSX.writetable("DOX_GEM.xlsx", df5)

    df6 = DataFrames.DataFrame(
        control_control = DOX_TAX[3, :, 1, 1],
        dox1taxol0_1 = DOX_TAX[3, :, 2, 2],
        dox1taxol1 = DOX_TAX[3, :, 2, 3],
        dox1taxol2 = DOX_TAX[3, :, 2, 4],
        dox1taxol3 = DOX_TAX[3, :, 2, 5],
        dox1taxol5 = DOX_TAX[3, :, 2, 6],
        dox1taxol7_5 = DOX_TAX[3, :, 2, 7],
        dox1taxol15 = DOX_TAX[3, :, 2, 8],
        dox10taxol0_1 = DOX_TAX[3, :, 3, 2],
        dox10taxol1 = DOX_TAX[3, :, 3, 3],
        dox10taxol2 = DOX_TAX[3, :, 3, 4],
        dox10taxol3 = DOX_TAX[3, :, 3, 5],
        dox10taxol5 = DOX_TAX[3, :, 3, 6],
        dox10taxol7_5 = DOX_TAX[3, :, 3, 7],
        dox10taxol15 = DOX_TAX[3, :, 3, 8],
        dox20taxol0_1 = DOX_TAX[3, :, 4, 2],
        dox20taxol1 = DOX_TAX[3, :, 4, 3],
        dox20taxol2 = DOX_TAX[3, :, 4, 4],
        dox20taxol3 = DOX_TAX[3, :, 4, 5],
        dox20taxol5 = DOX_TAX[3, :, 4, 6],
        dox20taxol7_5 = DOX_TAX[3, :, 4, 7],
        dox20taxol15 = DOX_TAX[3, :, 4, 8],
        dox50taxol0_1 = DOX_TAX[3, :, 5, 2],
        dox50taxol1 = DOX_TAX[3, :, 5, 3],
        dox50taxol2 = DOX_TAX[3, :, 5, 4],
        dox50taxol3 = DOX_TAX[3, :, 5, 5],
        dox50taxol5 = DOX_TAX[3, :, 5, 6],
        dox50taxol7_5 = DOX_TAX[3, :, 5, 7],
        dox50taxol15 = DOX_TAX[3, :, 5, 8],
        dox125taxol0_1 = DOX_TAX[3, :, 6, 2],
        dox125taxol1 = DOX_TAX[3, :, 6, 3],
        dox125taxol2 = DOX_TAX[3, :, 6, 4],
        dox125taxol3 = DOX_TAX[3, :, 6, 5],
        dox125taxol5 = DOX_TAX[3, :, 6, 6],
        dox125taxol7_5 = DOX_TAX[3, :, 6, 7],
        dox125taxol15 = DOX_TAX[3, :, 6, 8],
        dox250taxol0_1 = DOX_TAX[3, :, 7, 2],
        dox250taxol1 = DOX_TAX[3, :, 7, 3],
        dox250taxol2 = DOX_TAX[3, :, 7, 4],
        dox250taxol3 = DOX_TAX[3, :, 7, 5],
        dox250taxol5 = DOX_TAX[3, :, 7, 6],
        dox250taxol7_5 = DOX_TAX[3, :, 7, 7],
        dox250taxol15 = DOX_TAX[3, :, 7, 8],
        dox500taxol0_1 = DOX_TAX[3, :, 8, 2],
        dox500taxol1 = DOX_TAX[3, :, 8, 3],
        dox500taxol2 = DOX_TAX[3, :, 8, 4],
        dox500taxol3 = DOX_TAX[3, :, 8, 5],
        dox500taxol5 = DOX_TAX[3, :, 8, 6],
        dox500taxol7_5 = DOX_TAX[3, :, 8, 7],
        dox500taxol15 = DOX_TAX[3, :, 8, 8],
    )
    XLSX.writetable("DOX_TAX.xlsx", df6)

    df7 = DataFrames.DataFrame(
        control_control = DOX_PLB[3, :, 1, 1],
        dox1palbo5 = DOX_PLB[3, :, 2, 2],
        dox1palbo10 = DOX_PLB[3, :, 2, 3],
        dox1palbo25 = DOX_PLB[3, :, 2, 4],
        dox1palbo50 = DOX_PLB[3, :, 2, 5],
        dox1palbo100 = DOX_PLB[3, :, 2, 6],
        dox1palbo250 = DOX_PLB[3, :, 2, 7],
        dox1palbo500 = DOX_PLB[3, :, 2, 8],
        dox10palbo5 = DOX_PLB[3, :, 3, 2],
        dox10palbo10 = DOX_PLB[3, :, 3, 3],
        dox10palbo25 = DOX_PLB[3, :, 3, 4],
        dox10palbo50 = DOX_PLB[3, :, 3, 5],
        dox10palbo100 = DOX_PLB[3, :, 3, 6],
        dox10palbo250 = DOX_PLB[3, :, 3, 7],
        dox10palbo500 = DOX_PLB[3, :, 3, 8],
        dox20palbo5 = DOX_PLB[3, :, 4, 2],
        dox20palbo10 = DOX_PLB[3, :, 4, 3],
        dox20palbo25 = DOX_PLB[3, :, 4, 4],
        dox20palbo50 = DOX_PLB[3, :, 4, 5],
        dox20palbo100 = DOX_PLB[3, :, 4, 6],
        dox20palbo250 = DOX_PLB[3, :, 4, 7],
        dox20palbo500 = DOX_PLB[3, :, 4, 8],
        dox50palbo5 = DOX_PLB[3, :, 5, 2],
        dox50palbo10 = DOX_PLB[3, :, 5, 3],
        dox50palbo25 = DOX_PLB[3, :, 5, 4],
        dox50palbo50 = DOX_PLB[3, :, 5, 5],
        dox50palbo100 = DOX_PLB[3, :, 5, 6],
        dox50palbo250 = DOX_PLB[3, :, 5, 7],
        dox50palbo500 = DOX_PLB[3, :, 5, 8],
        dox125palbo5 = DOX_PLB[3, :, 6, 2],
        dox125palbo10 = DOX_PLB[3, :, 6, 3],
        dox125palbo25 = DOX_PLB[3, :, 6, 4],
        dox125palbo50 = DOX_PLB[3, :, 6, 5],
        dox125palbo100 = DOX_PLB[3, :, 6, 6],
        dox125palbo250 = DOX_PLB[3, :, 6, 7],
        dox125palbo500 = DOX_PLB[3, :, 6, 8],
        dox250palbo5 = DOX_PLB[3, :, 7, 2],
        dox250palbo10 = DOX_PLB[3, :, 7, 3],
        dox250palbo25 = DOX_PLB[3, :, 7, 4],
        dox250palbo50 = DOX_PLB[3, :, 7, 5],
        dox250palbo100 = DOX_PLB[3, :, 7, 6],
        dox250palbo250 = DOX_PLB[3, :, 7, 7],
        dox250palbo500 = DOX_PLB[3, :, 7, 8],
        dox500palbo5 = DOX_PLB[3, :, 8, 2],
        dox500palbo10 = DOX_PLB[3, :, 8, 3],
        dox500palbo25 = DOX_PLB[3, :, 8, 4],
        dox500palbo50 = DOX_PLB[3, :, 8, 5],
        dox500palbo100 = DOX_PLB[3, :, 8, 6],
        dox500palbo250 = DOX_PLB[3, :, 8, 7],
        dox500palbo500 = DOX_PLB[3, :, 8, 8],
    )

    XLSX.writetable("DOX_PLB.xlsx", df7)

    df8 = DataFrames.DataFrame(
        control_control = GEM_TAX[3, :, 1, 1],
        gem0_25taxol0_1 = GEM_TAX[3, :, 2, 2],
        gem0_25taxol1 = GEM_TAX[3, :, 2, 3],
        gem0_25taxol2 = GEM_TAX[3, :, 2, 4],
        gem0_25taxol3 = GEM_TAX[3, :, 2, 5],
        gem0_25taxol5 = GEM_TAX[3, :, 2, 6],
        gem0_25taxol7_5 = GEM_TAX[3, :, 2, 7],
        gem0_25taxol15 = GEM_TAX[3, :, 2, 8],
        gem1taxol0_1 = GEM_TAX[3, :, 3, 2],
        gem1taxol1 = GEM_TAX[3, :, 3, 3],
        gem1taxol2 = GEM_TAX[3, :, 3, 4],
        gem1taxol3 = GEM_TAX[3, :, 3, 5],
        gem1taxol5 = GEM_TAX[3, :, 3, 6],
        gem1taxol7_5 = GEM_TAX[3, :, 3, 7],
        gem1taxol15 = GEM_TAX[3, :, 3, 8],
        gem2_5taxol0_1 = GEM_TAX[3, :, 4, 2],
        gem2_5taxol1 = GEM_TAX[3, :, 4, 3],
        gem2_5taxol2 = GEM_TAX[3, :, 4, 4],
        gem2_5taxol3 = GEM_TAX[3, :, 4, 5],
        gem2_5taxol5 = GEM_TAX[3, :, 4, 6],
        gem2_5taxol7_5 = GEM_TAX[3, :, 4, 7],
        gem2_5taxol15 = GEM_TAX[3, :, 4, 8],
        gem5taxol0_1 = GEM_TAX[3, :, 5, 2],
        gem5taxol1 = GEM_TAX[3, :, 5, 3],
        gem5taxol2 = GEM_TAX[3, :, 5, 4],
        gem5taxol3 = GEM_TAX[3, :, 5, 5],
        gem5taxol5 = GEM_TAX[3, :, 5, 6],
        gem5taxol7_5 = GEM_TAX[3, :, 5, 7],
        gem5taxol15 = GEM_TAX[3, :, 5, 8],
        gem10taxol0_1 = GEM_TAX[3, :, 6, 2],
        gem10taxol1 = GEM_TAX[3, :, 6, 3],
        gem10taxol2 = GEM_TAX[3, :, 6, 4],
        gem10taxol3 = GEM_TAX[3, :, 6, 5],
        gem10taxol5 = GEM_TAX[3, :, 6, 6],
        gem10taxol7_5 = GEM_TAX[3, :, 6, 7],
        gem10taxol15 = GEM_TAX[3, :, 6, 8],
        gem17taxol0_1 = GEM_TAX[3, :, 7, 2],
        gem17taxol1 = GEM_TAX[3, :, 7, 3],
        gem17taxol2 = GEM_TAX[3, :, 7, 4],
        gem17taxol3 = GEM_TAX[3, :, 7, 5],
        gem17taxol5 = GEM_TAX[3, :, 7, 6],
        gem17taxol7_5 = GEM_TAX[3, :, 7, 7],
        gem17taxol15 = GEM_TAX[3, :, 7, 8],
        gem30taxol0_1 = GEM_TAX[3, :, 8, 2],
        gem30taxol1 = GEM_TAX[3, :, 8, 3],
        gem30taxol2 = GEM_TAX[3, :, 8, 4],
        gem30taxol3 = GEM_TAX[3, :, 8, 5],
        gem30taxol5 = GEM_TAX[3, :, 8, 6],
        gem30taxol7_5 = GEM_TAX[3, :, 8, 7],
        gem30taxol15 = GEM_TAX[3, :, 8, 8],
    )
    XLSX.writetable("GEM_TAX.xlsx", df8)

    df9 = DataFrames.DataFrame(
        control_control = GEM_PLB[3, :, 1, 1],
        gem0_25palbo5 = GEM_PLB[3, :, 2, 2],
        gem0_25palbo10 = GEM_PLB[3, :, 2, 3],
        gem0_25palbo25 = GEM_PLB[3, :, 2, 4],
        gem0_25palbo50 = GEM_PLB[3, :, 2, 5],
        gem0_25palbo100 = GEM_PLB[3, :, 2, 6],
        gem0_25palbo250 = GEM_PLB[3, :, 2, 7],
        gem0_25palbo500 = GEM_PLB[3, :, 2, 8],
        gem1palbo5 = GEM_PLB[3, :, 3, 2],
        gem1palbo10 = GEM_PLB[3, :, 3, 3],
        gem1palbo25 = GEM_PLB[3, :, 3, 4],
        gem1palbo50 = GEM_PLB[3, :, 3, 5],
        gem1palbo100 = GEM_PLB[3, :, 3, 6],
        gem1palbo250 = GEM_PLB[3, :, 3, 7],
        gem1palbo500 = GEM_PLB[3, :, 3, 8],
        gem2_5palbo5 = GEM_PLB[3, :, 4, 2],
        gem2_5palbo10 = GEM_PLB[3, :, 4, 3],
        gem2_5palbo25 = GEM_PLB[3, :, 4, 4],
        gem2_5palbo50 = GEM_PLB[3, :, 4, 5],
        gem2_5palbo100 = GEM_PLB[3, :, 4, 6],
        gem2_5palbo250 = GEM_PLB[3, :, 4, 7],
        gem2_5palbo500 = GEM_PLB[3, :, 4, 8],
        gem5palbo5 = GEM_PLB[3, :, 5, 2],
        gem5palbo10 = GEM_PLB[3, :, 5, 3],
        gem5palbo25 = GEM_PLB[3, :, 5, 4],
        gem5palbo50 = GEM_PLB[3, :, 5, 5],
        gem5palbo100 = GEM_PLB[3, :, 5, 6],
        gem5palbo250 = GEM_PLB[3, :, 5, 7],
        gem5palbo500 = GEM_PLB[3, :, 5, 8],
        gem10palbo5 = GEM_PLB[3, :, 6, 2],
        gem10palbo10 = GEM_PLB[3, :, 6, 3],
        gem10palbo25 = GEM_PLB[3, :, 6, 4],
        gem10palbo50 = GEM_PLB[3, :, 6, 5],
        gem10palbo100 = GEM_PLB[3, :, 6, 6],
        gem10palbo250 = GEM_PLB[3, :, 6, 7],
        gem10palbo500 = GEM_PLB[3, :, 6, 8],
        gem17palbo5 = GEM_PLB[3, :, 7, 2],
        gem17palbo10 = GEM_PLB[3, :, 7, 3],
        gem17palbo25 = GEM_PLB[3, :, 7, 4],
        gem17palbo50 = GEM_PLB[3, :, 7, 5],
        gem17palbo100 = GEM_PLB[3, :, 7, 6],
        gem17palbo250 = GEM_PLB[3, :, 7, 7],
        gem17palbo500 = GEM_PLB[3, :, 7, 8],
        gem30palbo5 = GEM_PLB[3, :, 8, 2],
        gem30palbo10 = GEM_PLB[3, :, 8, 3],
        gem30palbo25 = GEM_PLB[3, :, 8, 4],
        gem30palbo50 = GEM_PLB[3, :, 8, 5],
        gem30palbo100 = GEM_PLB[3, :, 8, 6],
        gem30palbo250 = GEM_PLB[3, :, 8, 7],
        gem30palbo500 = GEM_PLB[3, :, 8, 8],
    )

    XLSX.writetable("GEM_PLB.xlsx", df9)

    df10 = DataFrames.DataFrame(
        control_control = TAX_PLB[3, :, 1, 1],
        tax0_1palbo5 = TAX_PLB[3, :, 2, 2],
        tax0_1palbo10 = TAX_PLB[3, :, 2, 3],
        tax0_1palbo25 = TAX_PLB[3, :, 2, 4],
        tax0_1palbo50 = TAX_PLB[3, :, 2, 5],
        tax0_1palbo100 = TAX_PLB[3, :, 2, 6],
        tax0_1palbo250 = TAX_PLB[3, :, 2, 7],
        tax0_1palbo500 = TAX_PLB[3, :, 2, 8],
        tax1palbo5 = TAX_PLB[3, :, 3, 2],
        tax1palbo10 = TAX_PLB[3, :, 3, 3],
        tax1palbo25 = TAX_PLB[3, :, 3, 4],
        tax1palbo50 = TAX_PLB[3, :, 3, 5],
        tax1palbo100 = TAX_PLB[3, :, 3, 6],
        tax1palbo250 = TAX_PLB[3, :, 3, 7],
        tax1palbo500 = TAX_PLB[3, :, 3, 8],
        tax2palbo5 = TAX_PLB[3, :, 4, 2],
        tax2palbo10 = TAX_PLB[3, :, 4, 3],
        tax2palbo25 = TAX_PLB[3, :, 4, 4],
        tax2palbo50 = TAX_PLB[3, :, 4, 5],
        tax2palbo100 = TAX_PLB[3, :, 4, 6],
        tax2palbo250 = TAX_PLB[3, :, 4, 7],
        tax2palbo500 = TAX_PLB[3, :, 4, 8],
        tax3palbo5 = TAX_PLB[3, :, 5, 2],
        tax3palbo10 = TAX_PLB[3, :, 5, 3],
        tax3palbo25 = TAX_PLB[3, :, 5, 4],
        tax3palbo50 = TAX_PLB[3, :, 5, 5],
        tax3palbo100 = TAX_PLB[3, :, 5, 6],
        tax3palbo250 = TAX_PLB[3, :, 5, 7],
        tax3palbo500 = TAX_PLB[3, :, 5, 8],
        tax5palbo5 = TAX_PLB[3, :, 6, 2],
        tax5palbo10 = TAX_PLB[3, :, 6, 3],
        tax5palbo25 = TAX_PLB[3, :, 6, 4],
        tax5palbo50 = TAX_PLB[3, :, 6, 5],
        tax5palbo100 = TAX_PLB[3, :, 6, 6],
        tax5palbo250 = TAX_PLB[3, :, 6, 7],
        tax5palbo500 = TAX_PLB[3, :, 6, 8],
        taxol7_5palbo5 = TAX_PLB[3, :, 7, 2],
        taxol7_5palbo10 = TAX_PLB[3, :, 7, 3],
        taxol7_5palbo25 = TAX_PLB[3, :, 7, 4],
        taxol7_5palbo50 = TAX_PLB[3, :, 7, 5],
        taxol7_5palbo100 = TAX_PLB[3, :, 7, 6],
        taxol7_5palbo250 = TAX_PLB[3, :, 7, 7],
        taxol7_5palbo500 = TAX_PLB[3, :, 7, 8],
        taxol15palbo5 = TAX_PLB[3, :, 8, 2],
        taxol15palbo10 = TAX_PLB[3, :, 8, 3],
        taxol15palbo25 = TAX_PLB[3, :, 8, 4],
        taxol15palbo50 = TAX_PLB[3, :, 8, 5],
        taxol15palbo100 = TAX_PLB[3, :, 8, 6],
        taxol15palbo250 = TAX_PLB[3, :, 8, 7],
        taxol15palbo500 = TAX_PLB[3, :, 8, 8],
    )

    XLSX.writetable("TAX_PLB.xlsx", df10)
end


function combined_phaseDurations()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    t = LinRange(0.0, 95.0, 189)
    p = parameters()
    efcs = getODEparams(p, concs)
    # gem17 = DrugResponseModel.find_gem17(p)
    # efcs[:, 8, 3] = efcs[:, 7, 3]
    # efcs[:, 7, 3] = gem17

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

    function easier_(efcs)
        gi = zeros(2, 8, 8)
        gi[1, :, :] .= (4 ./ efcs[1, :, :] .+ 4 ./ efcs[2, :, :])
        gi[2, :, :] .= (5 ./ efcs[3, :, :] .+ 5 ./ efcs[4, :, :] .+ 5 ./ efcs[5, :, :] .+ 5 ./ efcs[6, :, :])
        gi
    end
    gim = zeros(2, 8, 8, 10)
    gim[:, :, :, 1] = easier_(LPT_DOX)
    gim[:, :, :, 2] = easier_(LPT_GEM)
    gim[:, :, :, 3] = easier_(LPT_TAX)
    gim[:, :, :, 4] = easier_(LPT_PLB)
    gim[:, :, :, 5] = easier_(DOX_GEM)
    gim[:, :, :, 6] = easier_(DOX_TAX)
    gim[:, :, :, 7] = easier_(DOX_PLB)
    gim[:, :, :, 8] = easier_(GEM_TAX)
    gim[:, :, :, 9] = easier_(GEM_PLB)
    gim[:, :, :, 10] = easier_(TAX_PLB)
    dfG1 = DataFrames.DataFrame(palbo50_lpt25=gim[1, 4, 5, 4], palbo50_lpt50=gim[1, 5, 5, 4], palbo50_lpt100=gim[1, 6, 5, 4], palbo50_lpt250=gim[1, 7, 5, 4], palbo50_gem5=gim[1, 5, 5, 9], palbo50_gem10=gim[1, 6, 5, 9], palbo50_gem30=gim[1, 7, 5, 9], palbo50_gem100=gim[1, 8, 5, 9], gem10_lpt25 = gim[1, 4, 6, 2], gem10_lpt50 = gim[1, 5, 6, 2], gem10_lpt100 = gim[1, 6, 6, 2], gem10_lpt250 = gim[1, 7, 6, 2], gem10_palbo25=gim[1, 6, 4, 9], gem10_palbo50=gim[1, 6, 5, 9], gem10_palbo100=gim[1, 6, 6, 9], gem10_palbo250=gim[1, 6, 7, 9], lap100_palbo25=gim[1, 6, 4, 4], lap100_palbo50=gim[1, 6, 5, 4], lap100_palbo100=gim[1, 6, 6, 4], lap100_palbo250=gim[1, 6, 7, 4], lap100_gem5=gim[1, 6, 5, 2], lap100_gem10=gim[1, 6, 6, 2], lap100_gem30=gim[1, 6, 7, 2], lap100_gem100=gim[1, 6, 8, 2], pax2_lpt25=gim[1, 4, 4, 3], pax2_lpt50=gim[1, 5, 4, 3], pax2_lpt100=gim[1, 6, 4, 3], pax2_lpt250=gim[1, 7, 4, 3])
    dfG2 = DataFrames.DataFrame(palbo50_lpt25=gim[2, 4, 5, 4], palbo50_lpt50=gim[2, 5, 5, 4], palbo50_lpt100=gim[2, 6, 5, 4], palbo50_lpt250=gim[2, 7, 5, 4], palbo50_gem5=gim[2, 5, 5, 9], palbo50_gem10=gim[2, 6, 5, 9], palbo50_gem30=gim[2, 7, 5, 9], palbo50_gem100=gim[2, 8, 5, 9], gem10_lpt25 = gim[2, 4, 6, 2], gem10_lpt50 = gim[2, 5, 6, 2], gem10_lpt100 = gim[2, 6, 6, 2], gem10_lpt250 = gim[2, 7, 6, 2], gem10_palbo25=gim[2, 6, 4, 9], gem10_palbo50=gim[2, 6, 5, 9], gem10_palbo100=gim[2, 6, 6, 9], gem10_palbo250=gim[2, 6, 7, 9], lap100_palbo25=gim[2, 6, 4, 4], lap100_palbo50=gim[2, 6, 5, 4], lap100_palbo100=gim[2, 6, 6, 4], lap100_palbo250=gim[2, 6, 7, 4], lap100_gem5=gim[2, 6, 5, 2], lap100_gem10=gim[2, 6, 6, 2], lap100_gem30=gim[2, 6, 7, 2], lap100_gem100=gim[2, 6, 8, 2], pax2_lpt25=gim[2, 4, 4, 3], pax2_lpt50=gim[2, 5, 4, 3], pax2_lpt100=gim[2, 6, 4, 3], pax2_lpt250=gim[2, 7, 4, 3])

    XLSX.writetable("G1CombinationDurations.xlsx", dfG1)
    XLSX.writetable("G2CombinationDurations.xlsx", dfG2)
end
