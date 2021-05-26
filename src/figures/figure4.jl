""" Figure 4: drug combination. """


function SSEs_combination()
    cellnum = DrugResponseModel.output_Bliss_cellnum() # bliss on cell number [189, 5, 5, 10]
    g = JLD.load("g.jld")["g"] # experimental data
    gT = g[1, 1:189, :, :] .+ g[2, 1:189, :, :] # [189, 6, 8]

    p = [0.151871, 0.677644, 1.24172, 1.22174, 0.0796602, 1.65385, 0.000103771, 0.000137389, 0.00010086, 0.000107739, 0.000100047, 0.0496591, 0.243334, 0.61755, 0.404729, 0.566801, 1.04495, 1.81665, 0.000103482, 0.000114637, 0.000101504, 0.000185468, 0.000101091, 0.199953, 63.3619, 2.91267, 0.000296213, 0.0350896, 1.37231, 0.0599936, 2.22239, 0.683019, 0.00978513, 0.000100827, 0.0419537, 0.000107985, 0.000116839, 0.0347002, 8.62214, 19.5609, 3.99883, 0.187608, 0.431694, 0.0710716, 3.17961, 0.550232, 0.00308626, 0.000101195, 0.000102833, 0.0066353, 0.0241068, 0.000104833, 22.9465, 0.621644, 3.05299, 0.156228, 0.361387, 1.03003, 0.229558, 0.807047, 0.000492391, 0.000100816, 0.000102552, 0.0165746, 0.000136098, 0.0587433, 3.9392, 0.228883, 1.10867, 0.24115, 2.52559, 0.288521]
    combinEffects = DrugResponseModel.my_helper(p)

    Gsim = zeros(2, 189, 6, 8)
    t = LinRange(0.0, 95.0, 189)

    for j=1:8
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
    for i=1:7
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

    p = [0.151871, 0.677644, 1.24172, 1.22174, 0.0796602, 1.65385, 0.000103771, 0.000137389, 0.00010086, 0.000107739, 0.000100047, 0.0496591, 0.243334, 0.61755, 0.404729, 0.566801, 1.04495, 1.81665, 0.000103482, 0.000114637, 0.000101504, 0.000185468, 0.000101091, 0.199953, 63.3619, 2.91267, 0.000296213, 0.0350896, 1.37231, 0.0599936, 2.22239, 0.683019, 0.00978513, 0.000100827, 0.0419537, 0.000107985, 0.000116839, 0.0347002, 8.62214, 19.5609, 3.99883, 0.187608, 0.431694, 0.0710716, 3.17961, 0.550232, 0.00308626, 0.000101195, 0.000102833, 0.0066353, 0.0241068, 0.000104833, 22.9465, 0.621644, 3.05299, 0.156228, 0.361387, 1.03003, 0.229558, 0.807047, 0.000492391, 0.000100816, 0.000102552, 0.0165746, 0.000136098, 0.0587433, 3.9392, 0.228883, 1.10867, 0.24115, 2.52559, 0.288521]
    combinEffects = DrugResponseModel.my_helper(p)

    Gsim = zeros(2, 189, 6, 8)
    t = LinRange(0.0, 95.0, 189)

    for j=1:8
        for i = 1:6 # concentration number
            Gsim[1, :, i, j], Gsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
        end
    end
    GsimT = Gsim[1, :, :, :] .+ Gsim[2, :, :, :] # Total model predictions    [189, 6, 8]

    SSEs = zeros(5, 2, 7) # dim1: Bliss on cell number - exp, dim2: Bliss on model - exp

    # avg(norm(Bliss on cell number - exp))
    SSEs[:, 1, 1] = sum((cellnum[1:189, :, 3, 4] .- gT[1:189, 2:6, 1]), dims=1) / 189
    SSEs[:, 1, 2] = sum((cellnum[1:189, :, 3, 9] .- gT[1:189, 2:6, 2]), dims=1) / 189
    SSEs[:, 1, 3] = sum((cellnum[1:189, 3, :, 9] .- gT[1:189, 2:6, 3]), dims=1) / 189
    SSEs[:, 1, 4] = sum((cellnum[1:189, :, 3, 2] .- gT[1:189, 2:6, 4]), dims=1) / 189
    SSEs[:, 1, 5] = sum((cellnum[1:189, 4, :, 4] .- gT[1:189, 2:6, 5]), dims=1) / 189
    SSEs[:, 1, 6] = sum((cellnum[1:189, 4, :, 2] .- gT[1:189, 2:6, 6]), dims=1) / 189
    SSEs[:, 1, 7] = sum((cellnum[1:189, 2, :, 5] .- gT[1:189, 2:6, 7]), dims=1) / 189

    # avg(norm(model predictions - exp))
    for i=1:7
        SSEs[:, 2, i] = sum((GsimT[1:189, 2:6, i] .- gT[1:189, 2:6, i]), dims=1) / 189
    end

    p1 = plot_sse_sep(repeat(["palbo 50 nM", "palbo 50 nM + lap25 nM", "palbo 50 nM + lap 50 nM", "palbo 50 nM + lap 100 nM", "palbo 50 nM + lap 250 nM"], outer=2), SSEs[:, :, 1])
    p2 = plot_sse_sep(repeat(["palbo 50 nM", "palbo 50 nM + gem 5 nM", "palbo 50 nM + gem 10 nM", "palbo 50 nM + gem 17 nM", "palbo 50 nM + gem 30 nM"], outer=2), SSEs[:, :, 2])
    p3 = plot_sse_sep(repeat(["gem 10 nM", "gem 10 nM + palbo 25 nM", "gem 10 nM + palbo 50 nM", "gem 10 nM + palbo 100 nM", "gem 10 nM + palbo 250 nM"], outer=2), SSEs[:, :, 3])
    p4 = plot_sse_sep(repeat(["gem 10 nM", "gem 10 nM + lap 25 nM", "gem 10 nM + lap 50 nM", "gem 10 nM + lap 100 nM", "gem 10 nM + lap 250 nM"], outer=2), SSEs[:, :, 4])
    p5 = plot_sse_sep(repeat(["lap 100 nM", "lap 100 nM + palbo 25 nM", "lap 100 nM + palbo 50 nM", "lap 100 nM + palbo 100 nM", "lap 100 nM + palbo 250 nM"], outer=2), SSEs[:, :, 5])
    p6 = plot_sse_sep(repeat(["lap 100 nM", "lap 100 nM + gem 5 nM", "lap 100 nM + gem 10 nM", "lap 100 nM + gem 17 nM", "lap 100 nM + gem 30 nM"], outer=2), SSEs[:, :, 6])
    p7 = plot_sse_sep(repeat(["dox 20 nM", "dox 20 nM + gem 5 nM", "dox 20 nM + gem 10 nM", "dox 20 nM + gem 17 nM", "dox 20 nM + gem 30 nM"], outer=2), SSEs[:, :, 7])
    p8 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, size = (2000, 700), layout=(2, 4))
    savefig(p, "fig4a.svg")
end
